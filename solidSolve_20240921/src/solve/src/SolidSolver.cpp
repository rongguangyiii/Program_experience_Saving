#include "grid/include/globalData.h"
#include "solve/include/SolidSolver.h"
#include <cmath>

SolidSolver::SolidSolver(const std::shared_ptr<GridBase> g, size_t nx, size_t ny) : residual_(0.0), grid_(g), caltime_(0),
u_(nx, std::vector<double>(ny, 0.0)), u_new_(nx, std::vector<double>(ny, 0.0)),
v_(nx, std::vector<double>(ny, 0.0)), v_new_(nx, std::vector<double>(ny, 0.0)),
urhs_(nx, std::vector<double>(ny, 0.0)), vrhs_(nx, std::vector<double>(ny, 0.0)),
Q_(nx, std::vector<double>(ny, 0.0)), R_(nx, std::vector<double>(ny, 0.0)),
S_(nx, std::vector<double>(ny, 0.0)), T_(nx, std::vector<double>(ny, 0.0)),
DFuDksi_(nx, std::vector<double>(ny, 0.0)), DFvDeta_(nx, std::vector<double>(ny, 0.0))
{
	double x_start_ = grid_->getXstart();
	double y_start_ = grid_->getYstart();
	double x_end_ = grid_->getXend();
	double y_end_ = grid_->getYend();
	nx_ = nx;
	ny_ = ny;
	x_step_ = (x_end_ - x_start_) / (nx - 1);
	y_step_ = (y_end_ - y_start_) / (ny - 1);
	E_ = 72000;//����ģ��: 72Gpa = 72000Mpa = 72000N/mm2
	miu_ = 0.33;//Possion Radiuo
	/*-------------------------------------------------------------
		Ҫ�����ȶ�������:dt<(dx*dx)/(2*D),D��ʾ��ɢϵ�����ȴ���ϵ����
		���������ȡ��������ѧ��P32 --��֥��ʽ(2-20)�е�E/(1-miu*miu)
	------------------------------------------------------------- */
	c1_ = E_ / (1 - miu_ * miu_);
	c2_ = (1 - miu_) / 2;
	c3_ = (1 + miu_) / 2;
	double c = 0.75;//����ϵ��������ʱ�䲽��С���ȶ���������
	dt_ = c * std::pow(x_step_, 2) / (2 * c1_);
};

void SolidSolver::solve()
{
	GlobalData::createFolder("../../result");
	//initialize();
	initialize_ref();
	timeAdvance();
}

void SolidSolver::initialize() {
	for (size_t i = 0; i < nx_; i++)
	{
		for (size_t j = 0; j < ny_; j++)
		{
			u_[i][j] = 0.0;
			v_[i][j] = 0.0;
		}
	}
	//�±߽�v=0�̶�;
	for (int i = 0; i < nx_; ++i) {
		v_[i][0] = 0.0;
	}
	//��߽� u=0�̶�;
	for (int j = 0; j < ny_; ++j) {
		u_[0][j] = 0.0;
	}
	//�ұ߽� u=1�̶�;
	for (int j = 0; j < ny_; ++j) {
		u_[nx_ - 1][j] = 1.0;
	}

	u_new_ = u_;
	v_new_ = v_;

}

void SolidSolver::initialize_ref() {
	const auto& Xcoord = grid_->getXCoord();
	for (size_t i = 0; i < nx_; i++)
	{
		for (size_t j = 0; j < ny_; j++)
		{
			u_[i][j] = Xcoord[i][j];
			v_[i][j] = 0.0;
		}
	}

	u_new_ = u_;
	v_new_ = v_;
}

void SolidSolver::timeAdvance()
{
	size_t timestep = 0;
	residual_ = 1.0;
	writeTecplotFile(timestep);
	while (residual_ > GlobalData::tolerance)
		//while (timestep < 300)
	{
		//if (caltime_ > GlobalData::ctrltime || timestep > GlobalData::ctrlstep)
		if (caltime_ > GlobalData::ctrltime || residual_ < GlobalData::tolerance || timestep > GlobalData::ctrlstep)
			break;
		++timestep;
		//1.ŷ��ʱ���ƽ�
		solveEulerTrans();
		//solveEuler();
		std::cout << "Timestep: " << timestep << " Residual: " << residual_ << std::endl;
		//2.��u����Ϊu_new��׼��������һ��
		u_ = u_new_;
		v_ = v_new_;
		writeTecplotFile(timestep);
	}
}

void SolidSolver::solveEuler()
{
	computeRHS();
	for (size_t i = 1; i < nx_ - 1; i++)
	{
		for (size_t j = 1; j < ny_ - 1; j++)
		{
			u_new_[i][j] = u_[i][j] + dt_ * c1_ * urhs_[i][j];
			v_new_[i][j] = v_[i][j] + dt_ * c1_ * vrhs_[i][j];
		}
	}
	extrapBoundary();
	computeResidual();
}
void SolidSolver::computeRHS()
{
	for (size_t i = 1; i < nx_ - 1; i++)
	{
		for (size_t j = 1; j < ny_ - 1; j++)
		{
			urhs_[i][j] = d2u_dx2(i, j) + c2_ * d2u_dy2(i, j) + c3_ * d2v_dxdy(i, j);
			vrhs_[i][j] = d2v_dy2(i, j) + c2_ * d2v_dx2(i, j) + c3_ * d2u_dxdy(i, j);
		}
	}

}

void SolidSolver::computeRHSTrans()
{
	for (size_t i = 1; i < nx_ - 1; i++)
	{
		for (size_t j = 1; j < ny_ - 1; j++)
		{
			urhs_[i][j] = dFu_dksi(i, j) + dFv_deta(i, j);
			vrhs_[i][j] = dGu_dksi(i, j) + dGv_deta(i, j);
		}
	}

}
void SolidSolver::solveEulerTrans()
{
	size_t index = 0;
	const auto& transVec = grid_->getCoefTransVec();
	computeRHSTrans();
	for (size_t i = 1; i < nx_ - 1; i++)
	{
		for (size_t j = 1; j < ny_ - 1; j++)
		{
			index = i * ny_ + j;
			u_new_[i][j] = u_[i][j] + dt_ * transVec[index].jacob() * c1_ * urhs_[i][j];
			v_new_[i][j] = v_[i][j] + dt_ * transVec[index].jacob() * c1_ * vrhs_[i][j];
		}
	}
	extrapBoundary();
	computeResidual();
}

void SolidSolver::computeResidual()
{
	residual_ = 0.0;
	for (int i = 1; i < nx_ - 1; ++i) {
		for (int j = 1; j < ny_ - 1; ++j) {
			residual_ += std::pow(urhs_[i][j], 2);
		}
	}
	residual_ = std::sqrt(residual_);
	if (residual_ > 1E10)
	{
		throw std::runtime_error("ʱ�䲽�����ó������㲻�ȶ���");
	}
}

void SolidSolver::extrapBoundary()
{
	/*	u��������,u�������Ѿ��̶���u_l=0.0,u_r=1.0*;
		v��������,v�����Ѿ��̶���v_d=0.0;*/
		//1 ��߽�
	for (int j = 0; j < ny_; ++j) {
		v_new_[0][j] = v_new_[1][j];
	}
	//1 �ϱ߽�
	for (int i = 0; i < nx_; ++i) {
		v_new_[i][ny_ - 1] = v_new_[i][ny_ - 2];
		u_new_[i][ny_ - 1] = u_new_[i][ny_ - 2];
	}
	//1 �±߽�
	for (int i = 0; i < nx_; ++i) {
		u_new_[i][0] = u_new_[i][1];
	}
	//1 �ұ߽�
	for (int j = 0; j < ny_; ++j) {
		v_new_[nx_ - 1][j] = v_new_[nx_ - 2][j];
	}
}

void SolidSolver::writeTecplotFile(size_t timestep) const {
	//if (timestep % GlobalData::ctrlout != 0)
	//	return;
	const auto& x = grid_->getXCoord();
	const auto& y = grid_->getYCoord();
	std::ofstream outFile("../../result/solution_" + std::to_string(timestep) + ".dat");
	outFile << "Title=\"Solution Data\"\n";
	//outFile << "Variables=\"X\", \"Y\", \"U\", \"V\", \"U_err\"\n";
	outFile << "Variables=\"X\", \"Y\", \"U\", \"dFu/dksi\", \"dFv/deta\", \"du/dx\", \"U_err\"\n";
	outFile << "Zone I = " << nx_ << ", J = " << ny_ << ", F = POINT" << ", STRANDID=1" << ", SOLUTIONTIME=" << timestep << std::endl;
	for (size_t j = 0; j < ny_; ++j) {
		for (size_t i = 0; i < nx_; ++i) {
			double ref = 1.0;
			double err = std::fabs(x[i][j] - u_[i][j]);
			//outFile << x[i][j] << " " << y[i][j] << " " << u_[i][j] << " " << v_[i][j] << " " << err << "\n";
			outFile << x[i][j] << " " << y[i][j] << " " << u_[i][j] << " " << DFuDksi_[i][j] << " " << DFvDeta_[i][j] << " " << Q_[i][j] << " " << err << "\n";
		}
	}
	outFile.close();
	std::cout << "�� " << timestep << " �����������;" << std::endl;

}

//Directly discretized----------------------------------------------------------------------------------------------------
double SolidSolver::d2u_dx2(size_t i, size_t j)
{
	double result = (u_[i + 1][j] - 2 * u_[i][j] + u_[i - 1][j]) / std::pow(x_step_, 2);
	return result;
}

double SolidSolver::d2u_dy2(size_t i, size_t j)
{
	double result = (u_[i][j + 1] - 2 * u_[i][j] + u_[i][j - 1]) / std::pow(y_step_, 2);
	return result;
}

double SolidSolver::d2u_dxdy(size_t i, size_t j)
{
	double result = (u_[i + 1][j + 1] + u_[i - 1][j - 1] - u_[i - 1][j + 1] - u_[i + 1][j - 1]) / (4 * y_step_ * x_step_);
	return result;
}
//v----------------------------------------------------------------------------------------------------
double SolidSolver::d2v_dx2(size_t i, size_t j)
{
	double result = (v_[i + 1][j] - 2 * v_[i][j] + v_[i - 1][j]) / std::pow(x_step_, 2);
	return result;
}

double SolidSolver::d2v_dy2(size_t i, size_t j)
{
	double result = (v_[i][j + 1] - 2 * v_[i][j] + v_[i][j - 1]) / std::pow(y_step_, 2);
	return result;
}

double SolidSolver::d2v_dxdy(size_t i, size_t j)
{
	double result = (v_[i + 1][j + 1] + v_[i - 1][j - 1] - v_[i - 1][j + 1] - v_[i + 1][j - 1]) / (4 * y_step_ * x_step_);
	return result;
}
//one order ----------------------------------------------------------------------------------------------------
double SolidSolver::du_dksi(size_t i, size_t j) {
	double result = 0;
	if (i == 0)
		result = (u_[i + 1][j] - u_[i][j]);
	else if (i == nx_ - 1)
		result = (u_[i][j] - u_[i - 1][j]);
	else
		result = (u_[i + 1][j] - u_[i - 1][j]) / 2;
	return result;
}
double SolidSolver::du_deta(size_t i, size_t j) {
	double result = 0;
	if (j == 0)
		result = (u_[i][j + 1] - u_[i][j]);
	else if (j == ny_ - 1)
		result = (u_[i][j] - u_[i][j - 1]);
	else
		result = (u_[i][j + 1] - u_[i][j - 1]) / 2;
	return result;
}
double SolidSolver::dv_dksi(size_t i, size_t j) {
	double result = 0;
	if (i == 0)
		result = (v_[i + 1][j] - v_[i][j]);
	else if (i == nx_ - 1)
		result = (v_[i][j] - v_[i - 1][j]);
	else
		result = (v_[i + 1][j] - v_[i - 1][j]) / 2;
	return result;
}
double SolidSolver::dv_deta(size_t i, size_t j) {
	double result = 0;
	if (j == 0)
		result = (v_[i][j + 1] - v_[i][j]);
	else if (j == ny_ - 1)
		result = (v_[i][j] - v_[i][j - 1]);
	else
		result = (v_[i][j + 1] - v_[i][j - 1]) / 2;
	return result;
}
//Q,R,S,T ----------------------------------------------------------------------------------------------------
double SolidSolver::Q(size_t i, size_t j, const CoefTrans& trans)//du_dx
{
	double result = trans.ksi_x() * du_dksi(i, j) + trans.eta_x() * du_deta(i, j);
	//result = 1.0;
	Q_[i][j] = result;
	Q_[0][0] = 1.0; Q_[nx_ - 1][ny_ - 1] = 1.0; Q_[0][ny_ - 1] = 1.0; Q_[nx_ - 1][0] = 1.0;
	return result;
}
double SolidSolver::R(size_t i, size_t j, const CoefTrans& trans)//du_dy
{
	double result = trans.ksi_y() * du_dksi(i, j) + trans.eta_y() * du_deta(i, j);
	//result = 0.0;
	R_[i][j] = result;
	return result;
}
double SolidSolver::S(size_t i, size_t j, const CoefTrans& trans)//dv_dx
{
	double result = trans.ksi_x() * dv_dksi(i, j) + trans.eta_x() * dv_deta(i, j);
	//result = 0.0;
	S_[i][j] = result;
	return result;
}
double SolidSolver::T(size_t i, size_t j, const CoefTrans& trans)//dv_dy
{
	double result = trans.ksi_y() * dv_dksi(i, j) + trans.eta_y() * dv_deta(i, j);
	//result = 0.0;
	T_[i][j] = result;
	return result;
}
//-------------------------------------------------------------------------------------------------------
double SolidSolver::Fu(size_t i, size_t j, const CoefTrans& trans)
{
	double ksix = trans.ksi_x();
	double ksiy = trans.ksi_y();
	double J = trans.jacob();
	double result = ksix / J * (Q(i, j, trans) + c3_ * T(i, j, trans)) + ksiy / J * (c2_ * R(i, j, trans));
	return result;
}
double SolidSolver::Fv(size_t i, size_t j, const CoefTrans& trans)
{
	double etax = trans.eta_x();
	double etay = trans.eta_y();
	double J = trans.jacob();
	double result = etax / J * (Q(i, j, trans) + c3_ * T(i, j, trans)) + etay / J * (c2_ * R(i, j, trans));
	return result;
}
double SolidSolver::Gu(size_t i, size_t j, const CoefTrans& trans)
{
	double ksix = trans.ksi_x();
	double ksiy = trans.ksi_y();
	double J = trans.jacob();
	double result = ksiy / J * T(i, j, trans) + ksix / J * (c2_ * S(i, j, trans) + c3_ * R(i, j, trans));
	return result;
}
double SolidSolver::Gv(size_t i, size_t j, const CoefTrans& trans)
{
	double etax = trans.eta_x();
	double etay = trans.eta_y();
	double J = trans.jacob();
	double result = etay / J * T(i, j, trans) + etax / J * (c2_ * S(i, j, trans) + c3_ * R(i, j, trans));
	return result;
}
//------------------------------------------------------------------------------------
double SolidSolver::dFu_dksi(size_t i, size_t j)
{
	size_t index = i * ny_ + j;
	const auto& trans = grid_->getCoefTransVec()[index];
	double result = (Fu(i + 1, j, trans) - Fu(i - 1, j, trans)) / 2;//������Ҫ������ʲô����㣿���û��ͨ�����ѵĻ���
	DFuDksi_[i][j] = result;
	return result;
}
double SolidSolver::dFv_deta(size_t i, size_t j)
{
	size_t index = i * ny_ + j;
	const auto& trans = grid_->getCoefTransVec()[index];
	double result = (Fv(i, j + 1, trans) - Fv(i, j - 1, trans)) / 2;
	DFvDeta_[i][j] = result;
	return result;
}
double SolidSolver::dGu_dksi(size_t i, size_t j)
{
	size_t index = i * ny_ + j;
	const auto& trans = grid_->getCoefTransVec()[index];
	double result = (Gu(i + 1, j, trans) - Gu(i - 1, j, trans)) / 2;
	return result;
}
double SolidSolver::dGv_deta(size_t i, size_t j)
{
	size_t index = i * ny_ + j;
	const auto& trans = grid_->getCoefTransVec()[index];
	double result = (Gv(i, j + 1, trans) - Gv(i, j - 1, trans)) / 2;
	return result;
}