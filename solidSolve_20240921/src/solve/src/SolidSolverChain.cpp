#include "grid/include/globalData.h"
#include "solve/include/SolidSolverChain.h"
#include <cmath>

SolidSolverChain::SolidSolverChain(const std::shared_ptr<GridBase> g, size_t nx, size_t ny) : residual_(0.0), grid_(g), caltime_(0),
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
	useRef_ = false;
};

void SolidSolverChain::solve()
{
	GlobalData::createFolder("../../resultChain");
	if (useRef_) {
		initialize_ref();//ʹ������ֵ��ʼ�����ڱ���,�߽粻������
	}
	else {
		initialize();
	}
	timeAdvance();
}

void SolidSolverChain::initialize() {
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

void SolidSolverChain::initialize_ref() {
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

void SolidSolverChain::timeAdvance()
{
	size_t timestep = 0;
	residual_ = 1.0;
	writeTecplotFile(timestep);
	//while (residual_ > GlobalData::tolerance)
	while (timestep < 300)
	{
		if (caltime_ > GlobalData::ctrltime || timestep > GlobalData::ctrlstep)
			//if (caltime_ > GlobalData::ctrltime || residual_ < GlobalData::tolerance || timestep > GlobalData::ctrlstep)
			break;
		++timestep;
		//1.ŷ��ʱ���ƽ�
		solveEulerTrans();
		std::cout << "Timestep: " << timestep << " Residual: " << residual_ << std::endl;
		//2.��u����Ϊu_new��׼��������һ��
		u_ = u_new_;
		v_ = v_new_;
		writeTecplotFile(timestep);
	}
}

void SolidSolverChain::computeRHSTrans()
{
	calPartialDer();
	//calPartialDer_ref();
	for (size_t i = 1; i < nx_ - 1; i++)
	{
		for (size_t j = 1; j < ny_ - 1; j++)
		{
			urhs_[i][j] = dFu_dksi(i, j) + dFv_deta(i, j);
			vrhs_[i][j] = dGu_dksi(i, j) + dGv_deta(i, j);
		}
	}

}
void SolidSolverChain::solveEulerTrans()
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

	if (useRef_) {
		/* ʹ������ֵ��ʼ�����ڱ���,�߽粻������ */
	}
	else {
		extrapBoundary();
	}
	computeResidual();
}

void SolidSolverChain::computeResidual()
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

void SolidSolverChain::extrapBoundary()
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

void SolidSolverChain::writeTecplotFile(size_t timestep) const {
	if (timestep % GlobalData::ctrlout != 0)
		return;
	const auto& x = grid_->getXCoord();
	const auto& y = grid_->getYCoord();
	std::ofstream outFile("../../resultHlaf/solution_" + std::to_string(timestep) + ".dat");
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
//one order ----------------------------------------------------------------------------------------------------
double SolidSolverChain::du_dksi(size_t i, size_t j) {
	double result = 0;
	if (i == 0)
		result = (u_[i + 1][j] - u_[i][j]);
	else if (i == nx_ - 1)
		result = (u_[i][j] - u_[i - 1][j]);
	else
		result = (u_[i + 1][j] - u_[i - 1][j]) / 2;
	return result;
}
double SolidSolverChain::du_deta(size_t i, size_t j) {
	double result = 0;
	if (j == 0)
		result = (u_[i][j + 1] - u_[i][j]);
	else if (j == ny_ - 1)
		result = (u_[i][j] - u_[i][j - 1]);
	else
		result = (u_[i][j + 1] - u_[i][j - 1]) / 2;
	return result;
}
double SolidSolverChain::dv_dksi(size_t i, size_t j) {
	double result = 0;
	if (i == 0)
		result = (v_[i + 1][j] - v_[i][j]);
	else if (i == nx_ - 1)
		result = (v_[i][j] - v_[i - 1][j]);
	else
		result = (v_[i + 1][j] - v_[i - 1][j]) / 2;
	return result;
}
double SolidSolverChain::dv_deta(size_t i, size_t j) {
	double result = 0;
	if (j == 0)
		result = (v_[i][j + 1] - v_[i][j]);
	else if (j == ny_ - 1)
		result = (v_[i][j] - v_[i][j - 1]);
	else
		result = (v_[i][j + 1] - v_[i][j - 1]) / 2;
	return result;
}
//Fu,Fv,Gu,Gv-------------------------------------------------------------------------------------------------------
double SolidSolverChain::Fu(size_t i, size_t j, const CoefTrans& trans)
{
	double ksix = trans.ksi_x();
	double ksiy = trans.ksi_y();
	double J = trans.jacob();
	//std::cout << "DEBUG: FU�����\t[" << i << ", " << j << "]ʹ�õĶ���ϵ��Ϊ:   " << ksix / J << " \t" << ksiy / J << std::endl;
	double result = ksix / J * (Q_[i][j] + c3_ * T_[i][j]) + ksiy / J * (c2_ * R_[i][j]);
	return result;
}
double SolidSolverChain::Fv(size_t i, size_t j, const CoefTrans& trans)
{
	double etax = trans.eta_x();
	double etay = trans.eta_y();
	double J = trans.jacob();
	double result = etax / J * (Q_[i][j] + c3_ * T_[i][j]) + etay / J * (c2_ * R_[i][j]);
	return result;
}
double SolidSolverChain::Gu(size_t i, size_t j, const CoefTrans& trans)
{
	double ksix = trans.ksi_x();
	double ksiy = trans.ksi_y();
	double J = trans.jacob();
	double result = ksiy / J * T_[i][j] + ksix / J * (c2_ * S_[i][j] + c3_ * R_[i][j]);
	return result;
}
double SolidSolverChain::Gv(size_t i, size_t j, const CoefTrans& trans)
{
	double etax = trans.eta_x();
	double etay = trans.eta_y();
	double J = trans.jacob();
	double result = etay / J * T_[i][j] + etax / J * (c2_ * S_[i][j] + c3_ * R_[i][j]);
	return result;
}
//Half Point value------------------------------------------------------------------------------------------
double SolidSolverChain::FuHalfPlus(size_t i, size_t j, const CoefTrans& trans) {
	double halfresult = (Fu(i, j, trans) + Fu(i + 1, j, trans)) / 2;
	return halfresult;
}
double SolidSolverChain::FuHalfMius(size_t i, size_t j, const CoefTrans& trans) {
	double halfresult = (Fu(i, j, trans) + Fu(i - 1, j, trans)) / 2;
	return halfresult;
}
double SolidSolverChain::FvHalfPlus(size_t i, size_t j, const CoefTrans& trans) {
	double halfresult = (Fv(i, j, trans) + Fv(i, j + 1, trans)) / 2;
	return halfresult;
}
double SolidSolverChain::FvHalfMius(size_t i, size_t j, const CoefTrans& trans) {
	double halfresult = (Fv(i, j, trans) + Fv(i, j - 1, trans)) / 2;
	return halfresult;
}
double SolidSolverChain::GuHalfPlus(size_t i, size_t j, const CoefTrans& trans) {
	double halfresult = (Gu(i, j, trans) + Gu(i + 1, j, trans)) / 2;
	return halfresult;
}
double SolidSolverChain::GuHalfMius(size_t i, size_t j, const CoefTrans& trans) {
	double halfresult = (Gu(i, j, trans) + Gu(i - 1, j, trans)) / 2;
	return halfresult;
}
double SolidSolverChain::GvHalfPlus(size_t i, size_t j, const CoefTrans& trans) {
	double halfresult = (Gv(i, j, trans) + Gv(i, j + 1, trans)) / 2;
	return halfresult;
}
double SolidSolverChain::GvHalfMius(size_t i, size_t j, const CoefTrans& trans) {
	double halfresult = (Gv(i, j, trans) + Gv(i, j - 1, trans)) / 2;
	return halfresult;
}
//dFu_dksi,dGu_dksi------------------------------------------------------------------------------------
double SolidSolverChain::dFu_dksi(size_t i, size_t j)
{
	size_t index = i * ny_ + j;
	const auto& trans = grid_->getCoefTransVec()[index];
	//std::cout << "DEBUG: ��ǰ�����Ϊ\t[" << i << ", " << j << "]:   " << trans.ksi_x() / trans.jacob() << " \t" << trans.ksi_y() / trans.jacob() << std::endl;
	double result = (FuHalfPlus(i, j, trans) - FuHalfMius(i, j, trans));
	DFuDksi_[i][j] = result;
	return result;
}
double SolidSolverChain::dFv_deta(size_t i, size_t j)
{
	size_t index = i * ny_ + j;
	const auto& trans = grid_->getCoefTransVec()[index];
	double result = (FvHalfPlus(i, j, trans) - FvHalfMius(i, j, trans));
	DFvDeta_[i][j] = result;
	return result;
}
double SolidSolverChain::dGu_dksi(size_t i, size_t j)
{
	size_t index = i * ny_ + j;
	const auto& trans = grid_->getCoefTransVec()[index];
	double result = (GuHalfPlus(i, j, trans) - GuHalfMius(i, j, trans));
	return result;
}
double SolidSolverChain::dGv_deta(size_t i, size_t j)
{
	size_t index = i * ny_ + j;
	const auto& trans = grid_->getCoefTransVec()[index];
	double result = (GvHalfPlus(i, j, trans) - GvHalfMius(i, j, trans));
	return result;
}
//Chain rule---------------------------------------------------------------------------
void  SolidSolverChain::calPartialDer()
{
	std::vector<std::string> dirname{ "u","v" };
	const auto& x = grid_->getXCoord();
	const auto& y = grid_->getYCoord();
	for (size_t i = 0; i < nx_; i++)
	{
		for (size_t j = 0; j < ny_; j++)
		{
			//du_dx,du_dy
			Q_[i][j] = Q(i, j);
			R_[i][j] = R(i, j);

			//dv_dx,dv_dy
			S_[i][j] = S(i, j);
			T_[i][j] = T(i, j);
		}
	}
}
void  SolidSolverChain::calPartialDer_ref()
{
	std::vector<std::string> dirname{ "u","v" };
	const auto& x = grid_->getXCoord();
	const auto& y = grid_->getYCoord();
	for (size_t i = 0; i < nx_; i++)
	{
		for (size_t j = 0; j < ny_; j++)
		{
			//du_dx,du_dy
			Q_[i][j] = 1.0;
			R_[i][j] = 0.0;

			//dv_dx,dv_dy
			S_[i][j] = 0.0;
			T_[i][j] = 0.0;
		}
	}
}
//Q,R,S,T ----------------------------------------------------------------------------------------------------
double SolidSolverChain::Q(size_t i, size_t j)//du_dx
{
	size_t index = i * ny_ + j;
	const auto& trans = grid_->getCoefTransVec()[index];
	//std::cout << "DEBUG: Q�����\t[" << i << ", " << j << "]ʹ�õĶ���ϵ��Ϊ:   " << trans.ksi_x() << " \t" << trans.eta_x() << std::endl;
	double result = trans.ksi_x() * du_dksi(i, j) + trans.eta_x() * du_deta(i, j);
	return result;
}
double SolidSolverChain::R(size_t i, size_t j)//du_dy
{
	size_t index = i * ny_ + j;
	const auto& trans = grid_->getCoefTransVec()[index];
	double result = trans.ksi_y() * du_dksi(i, j) + trans.eta_y() * du_deta(i, j);
	return result;
}
double SolidSolverChain::S(size_t i, size_t j)//dv_dx
{
	size_t index = i * ny_ + j;
	const auto& trans = grid_->getCoefTransVec()[index];
	double result = trans.ksi_x() * dv_dksi(i, j) + trans.eta_x() * dv_deta(i, j);
	return result;
}
double SolidSolverChain::T(size_t i, size_t j)//dv_dy
{
	size_t index = i * ny_ + j;
	const auto& trans = grid_->getCoefTransVec()[index];
	double result = trans.ksi_y() * dv_dksi(i, j) + trans.eta_y() * dv_deta(i, j);
	return result;
}