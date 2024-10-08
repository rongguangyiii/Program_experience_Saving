#include "grid/include/globalData.h"
#include "solve/include/SolidSolverLSGNoDEER.h"
#include <cmath>

SolidSolverLSGNoDEER::SolidSolverLSGNoDEER(const std::shared_ptr<GridBase> g, size_t nx, size_t ny) : residual_(0.0), grid_(g), caltime_(0),
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
	E_ = 72000;//弹性模量: 72Gpa = 72000Mpa = 72000N/mm2
	miu_ = 0.33;//Possion Radiuo
	/*-------------------------------------------------------------
		要满足稳定性条件:dt<(dx*dx)/(2*D),D表示扩散系数或热传导系数。
		在这里可以取《弹性力学》P32 --徐芝纶式(2-20)中的E/(1-miu*miu)
	------------------------------------------------------------- */
	c1_ = E_ / (1 - miu_ * miu_);
	c2_ = (1 - miu_) / 2;
	c3_ = (1 + miu_) / 2;
	double c = 0.75;//控制系数，控制时间步长小于稳定性条件。
	dt_ = c * std::pow(x_step_, 2) / (2 * c1_);
	useRef_ = true;
};

void SolidSolverLSGNoDEER::solve()
{
	GlobalData::createFolder("../../resultLSGNoDeer");
	if (useRef_) {
		initialize_ref();//使用理论值初始化场内变量,边界不用外推
	}
	else {
		initialize();
	}
	timeAdvance();
}

void SolidSolverLSGNoDEER::initialize() {
	for (size_t i = 0; i < nx_; i++)
	{
		for (size_t j = 0; j < ny_; j++)
		{
			u_[i][j] = 0.0;
			v_[i][j] = 0.0;
		}
	}
	//下边界v=0固定;
	for (int i = 0; i < nx_; ++i) {
		v_[i][0] = 0.0;
	}
	//左边界 u=0固定;
	for (int j = 0; j < ny_; ++j) {
		u_[0][j] = 0.0;
	}
	//右边界 u=1固定;
	for (int j = 0; j < ny_; ++j) {
		u_[nx_ - 1][j] = 0.1 + j / 21.0;
		//u_[nx_ - 1][j] = 1.0;
	}

	u_new_ = u_;
	v_new_ = v_;

}

void SolidSolverLSGNoDEER::initialize_ref() {
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

void SolidSolverLSGNoDEER::timeAdvance()
{
	size_t timestep = 0;
	residual_ = 1.0;
	writeTecplotFile(timestep);
	//while (residual_ > GlobalData::tolerance)
	while (timestep < 100)
	{
		//if (caltime_ > GlobalData::ctrltime || residual_ < GlobalData::tolerance || timestep > GlobalData::ctrlstep)
		if (caltime_ > GlobalData::ctrltime || timestep > GlobalData::ctrlstep)
			break;
		++timestep;
		//1.欧拉时间推进
		solveEulerTrans();
		std::cout << "Timestep: " << timestep << " Residual: " << residual_ << std::endl;
		//2.将u更新为u_new，准备进入下一步
		u_ = u_new_;
		v_ = v_new_;
		writeTecplotFile(timestep);
	}
}

void SolidSolverLSGNoDEER::computeRHSTrans()
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
void SolidSolverLSGNoDEER::solveEulerTrans()
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
		/* 使用理论值初始化场内变量,边界不用外推 */
	}
	else {
		extrapBoundary();
	}
	computeResidual();
}

void SolidSolverLSGNoDEER::computeResidual()
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
		throw std::runtime_error("时间步长设置出错，计算不稳定！");
	}
}

void SolidSolverLSGNoDEER::extrapBoundary()
{
	/*	u上下左右,u的左右已经固定：u_l=0.0,u_r=1.0*;
		v上下左右,v的下已经固定：v_d=0.0;*/
		//1 左边界
	for (int j = 0; j < ny_; ++j) {
		v_new_[0][j] = v_new_[1][j];
	}
	//1 上边界
	for (int i = 0; i < nx_; ++i) {
		v_new_[i][ny_ - 1] = v_new_[i][ny_ - 2];
		u_new_[i][ny_ - 1] = u_new_[i][ny_ - 2];
	}
	//1 下边界
	for (int i = 0; i < nx_; ++i) {
		u_new_[i][0] = u_new_[i][1];
	}
	//1 右边界
	for (int j = 0; j < ny_; ++j) {
		v_new_[nx_ - 1][j] = v_new_[nx_ - 2][j];
	}
}

void SolidSolverLSGNoDEER::writeTecplotFile(size_t timestep) const {
	//if (timestep % GlobalData::ctrlout != 0)
	//	return;
	const auto& x = grid_->getXCoord();
	const auto& y = grid_->getYCoord();
	std::ofstream outFile("../../resultLSGNoDeer/solution_" + std::to_string(timestep) + ".dat");
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
	std::cout << "第 " << timestep << " 步结果输出完成;" << std::endl;

}
//Fu,Fv,Gu,Gv-------------------------------------------------------------------------------------------------------
double SolidSolverLSGNoDEER::Fu(size_t i, size_t j)
{
	size_t index = i * ny_ + j;
	const auto& trans = grid_->getCoefTransVec()[index];
	double ksix = trans.ksi_x();
	double ksiy = trans.ksi_y();
	double J = trans.jacob();
	//std::cout << "DEBUG: FU计算点\t[" << i << ", " << j << "]使用的度量系数为:   " << ksix / J << " \t" << ksiy / J << std::endl;
	double result = ksix / J * (Q_[i][j] + c3_ * T_[i][j]) + ksiy / J * (c2_ * R_[i][j]);
	return result;
}
double SolidSolverLSGNoDEER::Fv(size_t i, size_t j)
{
	size_t index = i * ny_ + j;
	const auto& trans = grid_->getCoefTransVec()[index];
	double etax = trans.eta_x();
	double etay = trans.eta_y();
	double J = trans.jacob();
	double result = etax / J * (Q_[i][j] + c3_ * T_[i][j]) + etay / J * (c2_ * R_[i][j]);
	return result;
}
double SolidSolverLSGNoDEER::Gu(size_t i, size_t j)
{
	size_t index = i * ny_ + j;
	const auto& trans = grid_->getCoefTransVec()[index];
	double ksix = trans.ksi_x();
	double ksiy = trans.ksi_y();
	double J = trans.jacob();
	double result = ksiy / J * T_[i][j] + ksix / J * (c2_ * S_[i][j] + c3_ * R_[i][j]);
	return result;
}
double SolidSolverLSGNoDEER::Gv(size_t i, size_t j)
{
	size_t index = i * ny_ + j;
	const auto& trans = grid_->getCoefTransVec()[index];
	double etax = trans.eta_x();
	double etay = trans.eta_y();
	double J = trans.jacob();
	double result = etay / J * T_[i][j] + etax / J * (c2_ * S_[i][j] + c3_ * R_[i][j]);
	return result;
}
//dFu_dksi,dFv_dksi------------------------------------------------------------------------------------
double SolidSolverLSGNoDEER::dFu_dksi(size_t i, size_t j)
{
	size_t index = i * ny_ + j;
	const auto& trans = grid_->getCoefTransVec()[index];
	//std::cout << "DEBUG: 当前计算点为\t[" << i << ", " << j << "]:   " << trans.ksi_x() / trans.jacob() << " \t" << trans.ksi_y() / trans.jacob() << std::endl;
	double result = (Fu(i, j) - Fu(i - 1, j)) / 2;
	//std::cout << "DEBUG:========================================================================================= \t" << std::endl;
	DFuDksi_[i][j] = result;
	return result;
}
double SolidSolverLSGNoDEER::dFv_deta(size_t i, size_t j)
{
	size_t index = i * ny_ + j;
	const auto& trans = grid_->getCoefTransVec()[index];
	double result = (Fv(i, j) - Fv(i, j - 1)) / 2;
	DFvDeta_[i][j] = result;
	return result;
}
double SolidSolverLSGNoDEER::dGu_dksi(size_t i, size_t j)
{
	size_t index = i * ny_ + j;
	const auto& trans = grid_->getCoefTransVec()[index];
	double result = (Gu(i, j) - Gu(i - 1, j)) / 2;
	return result;
}
double SolidSolverLSGNoDEER::dGv_deta(size_t i, size_t j)
{
	size_t index = i * ny_ + j;
	const auto& trans = grid_->getCoefTransVec()[index];
	double result = (Gv(i, j) - Gv(i, j - 1)) / 2;
	return result;
}
//LSGrad---------------------------------------------------------------------------
void  SolidSolverLSGNoDEER::calPartialDer()
{
	std::vector<std::string> dirname{ "u","v" };
	const auto& x = grid_->getXCoord();
	const auto& y = grid_->getYCoord();
	for (size_t i = 0; i < nx_; i++)
	{
		for (size_t j = 0; j < ny_; j++)
		{
			//du_dx,du_dy
			std::vector<LSGradCoord> neighborsVec_1 = genLSGradCoord(i, j, "u");
			LSGradCoord curlsGCoord_1(x[i][j], y[i][j], u_[i][j]);
			LSGrad dudx_dudy(neighborsVec_1, curlsGCoord_1);
			double ux = 0;
			double uy = 0;
			dudx_dudy.computeGradient(ux, uy);
			Q_[i][j] = ux;
			R_[i][j] = uy;

			//dv_dx,dv_dy
			std::vector<LSGradCoord> neighborsVec_2 = genLSGradCoord(i, j, "v");
			LSGradCoord curlsGCoord_2(x[i][j], y[i][j], v_[i][j]);
			LSGrad dvdx_dvdy(neighborsVec_2, curlsGCoord_2);
			double vx = 0;
			double vy = 0;
			dvdx_dvdy.computeGradient(vx, vy);
			S_[i][j] = vx;
			T_[i][j] = vy;
		}
	}
}
void  SolidSolverLSGNoDEER::calPartialDer_ref()
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
std::vector<LSGradCoord>  SolidSolverLSGNoDEER::genLSGradCoord(size_t i, size_t j, std::string flag)
{
	const auto& x = grid_->getXCoord();
	const auto& y = grid_->getYCoord();
	std::vector<std::vector<double>> val;
	if (flag == "u")
		val = u_;
	else if (flag == "v")
		val = v_;
	else
		throw std::runtime_error("给定方向出错！");
	std::vector<LSGradCoord> neighborsVec;
	if (i == 0 && j != 0 && j != ny_ - 1) {  // 1. 左边界 (i == 0) 非角点 (0 < j < my)
		neighborsVec.emplace_back(x[i + 1][j], y[i + 1][j], val[i + 1][j]);
		//neighborsVec.emplace_back(x[i - 1][j], y[i - 1][j], val[i - 1][j]);
		neighborsVec.emplace_back(x[i][j + 1], y[i][j + 1], val[i][j + 1]);
		neighborsVec.emplace_back(x[i][j - 1], y[i][j - 1], val[i][j - 1]);
	}
	else if (i == nx_ - 1 && j != 0 && j != ny_ - 1) {  // 2. 右边界 (i == mx) 非角点 (0 < j < my)
		//neighborsVec.emplace_back(x[i + 1][j], y[i + 1][j], val[i + 1][j]);
		neighborsVec.emplace_back(x[i - 1][j], y[i - 1][j], val[i - 1][j]);
		neighborsVec.emplace_back(x[i][j + 1], y[i][j + 1], val[i][j + 1]);
		neighborsVec.emplace_back(x[i][j - 1], y[i][j - 1], val[i][j - 1]);
	}
	else if (i == 0 && j == 0) {  // 3. 左下角 (i == 0, j == 0)
		neighborsVec.emplace_back(x[i + 1][j], y[i + 1][j], val[i + 1][j]);
		//neighborsVec.emplace_back(x[i - 1][j], y[i - 1][j], val[i - 1][j]);
		neighborsVec.emplace_back(x[i][j + 1], y[i][j + 1], val[i][j + 1]);
		//neighborsVec.emplace_back(x[i][j - 1], y[i][j - 1], val[i][j - 1]);
	}
	else if (i == 0 && j == ny_ - 1) {  // 4. 左上角 (i == 0, j == my)
		neighborsVec.emplace_back(x[i + 1][j], y[i + 1][j], val[i + 1][j]);
		//neighborsVec.emplace_back(x[i - 1][j], y[i - 1][j], val[i - 1][j]);
		//neighborsVec.emplace_back(x[i][j + 1], y[i][j + 1], val[i][j + 1]);
		neighborsVec.emplace_back(x[i][j - 1], y[i][j - 1], val[i][j - 1]);
	}
	else if (i == nx_ - 1 && j == ny_ - 1) {  // 5. 右上角 (i == mx, j == my)
		//neighborsVec.emplace_back(x[i + 1][j], y[i + 1][j], val[i + 1][j]);
		neighborsVec.emplace_back(x[i - 1][j], y[i - 1][j], val[i - 1][j]);
		//neighborsVec.emplace_back(x[i][j + 1], y[i][j + 1], val[i][j + 1]);
		neighborsVec.emplace_back(x[i][j - 1], y[i][j - 1], val[i][j - 1]);
	}
	else if (i == nx_ - 1 && j == 0) {  // 6. 右下角 (i == mx, j == 0)
		//neighborsVec.emplace_back(x[i + 1][j], y[i + 1][j], val[i + 1][j]);
		neighborsVec.emplace_back(x[i - 1][j], y[i - 1][j], val[i - 1][j]);
		neighborsVec.emplace_back(x[i][j + 1], y[i][j + 1], val[i][j + 1]);
		//neighborsVec.emplace_back(x[i][j - 1], y[i][j - 1], val[i][j - 1]);
	}
	else if (j == 0) {  // 7. 下边界 (j == 0) 非角点 (0 < i < mx)
		neighborsVec.emplace_back(x[i + 1][j], y[i + 1][j], val[i + 1][j]);
		neighborsVec.emplace_back(x[i - 1][j], y[i - 1][j], val[i - 1][j]);
		neighborsVec.emplace_back(x[i][j + 1], y[i][j + 1], val[i][j + 1]);
		//neighborsVec.emplace_back(x[i][j - 1], y[i][j - 1], val[i][j - 1]);
	}
	else if (j == ny_ - 1) {  // 8. 上边界 (j == my) 非角点 (0 < i < mx)
		neighborsVec.emplace_back(x[i + 1][j], y[i + 1][j], val[i + 1][j]);
		neighborsVec.emplace_back(x[i - 1][j], y[i - 1][j], val[i - 1][j]);
		//neighborsVec.emplace_back(x[i][j + 1], y[i][j + 1], val[i][j + 1]);
		neighborsVec.emplace_back(x[i][j - 1], y[i][j - 1], val[i][j - 1]);
	}
	else if (i != 0 && j != 0 && i != nx_ - 1 && j != ny_ - 1) {  // 9. 内部点
		neighborsVec.emplace_back(x[i + 1][j], y[i + 1][j], val[i + 1][j]);
		neighborsVec.emplace_back(x[i - 1][j], y[i - 1][j], val[i - 1][j]);
		neighborsVec.emplace_back(x[i][j + 1], y[i][j + 1], val[i][j + 1]);
		neighborsVec.emplace_back(x[i][j - 1], y[i][j - 1], val[i][j - 1]);
	}
	else
	{
		throw std::runtime_error("当前点的邻居出错！");
	}

	return neighborsVec;
}