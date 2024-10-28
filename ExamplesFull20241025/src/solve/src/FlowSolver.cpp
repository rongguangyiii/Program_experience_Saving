#include "grid/include/globalData.h"
#include "solve/include/flowSolver.h"
#include "tools/include/tools.h"
#include <cmath>

FlowSolver::FlowSolver(const std::shared_ptr<GridBase> g, size_t nx, size_t ny) : residual_(0.0), grid_(g), caltime_(0),
u_(nx, std::vector<double>(ny, 0.0)), u_new_(nx, std::vector<double>(ny, 0.0)),
urhs_(nx, std::vector<double>(ny, 0.0)), Q_(nx, std::vector<double>(ny, 0.0)), R_(nx, std::vector<double>(ny, 0.0)),
S_(nx, std::vector<double>(ny, 0.0)), T_(nx, std::vector<double>(ny, 0.0)),
dFdksi_(nx, std::vector<double>(ny, 0.0)), dGdeta_(nx, std::vector<double>(ny, 0.0))
{
	folder_name_ = "result";
	double x_start_ = grid_->getXstart();
	double y_start_ = grid_->getYstart();
	double x_end_ = grid_->getXend();
	double y_end_ = grid_->getYend();
	nx_ = nx;
	ny_ = ny;
	x_step_ = (x_end_ - x_start_) / (nx - 1);
	y_step_ = (y_end_ - y_start_) / (ny - 1);
	double c1 =1.0;
	double c = 0.7;//控制系数，控制时间步长小于稳定性条件。
	dt_ = c * x_step_ / c1;
	useRef_ = true;
};

void FlowSolver::solve()
{
	Tools::createFolder("../../" + folder_name_);
	if (useRef_) {
		initialize_ref();//使用理论值初始化场内变量,边界不用外推
	}
	else {
		initialize();
	}
	timeAdvance();
}

void FlowSolver::initialize() {
	for (size_t i = 0; i < nx_; i++)
	{
		for (size_t j = 0; j < ny_; j++)
		{
			u_[i][j] = 0.0;
		}
	}
	//左边界 u=0固定;
	for (int j = 0; j < ny_; ++j) {
		u_[0][j] = 1.0;
	}

	u_new_ = u_;
}

void FlowSolver::initialize_ref() {
	//const auto& Xcoord = grid_->getXCoord();
	for (size_t i = 0; i < nx_; i++)
	{
		for (size_t j = 0; j < ny_; j++)
		{
			u_[i][j] = 1.0;
		}
	}
	u_new_ = u_;
}

void FlowSolver::timeAdvance()
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
		//1.欧拉时间推进
		solveEulerTrans();
		//solveEuler();
		std::cout << "Timestep: " << timestep << " Residual: " << residual_ << std::endl;
		//2.将u更新为u_new，准备进入下一步
		u_ = u_new_;
		writeTecplotFile(timestep);
	}
}

void FlowSolver::solveEuler()
{
	computeRHS();
	for (size_t i = 1; i < nx_ - 1; i++)
	{
		for (size_t j = 1; j < ny_ - 1; j++)
		{
			u_new_[i][j] = u_[i][j] - dt_ * urhs_[i][j];
		}
	}
	extrapBoundary();
	computeResidual();
}

void FlowSolver::computeRHS()
{
	for (size_t i = 1; i < nx_ - 1; i++)
	{
		for (size_t j = 1; j < ny_ - 1; j++)
		{
			urhs_[i][j] = dau_dx(i, j) + dbu_dy(i, j);
		}
	}

}

void FlowSolver::solveEulerTrans()
{
	size_t index = 0;
	const auto& transVec = grid_->getCoefTransVec();
	computeRHSTrans();
	for (size_t i = 1; i < nx_ - 1; i++)
	{
		for (size_t j = 1; j < ny_ - 1; j++)
		{
			index = i * ny_ + j;
			u_new_[i][j] = u_[i][j] - dt_ * transVec[index].jacob() * urhs_[i][j];
		}
	}
	extrapBoundary();
	computeResidual();
}

void FlowSolver::computeRHSTrans()
{
	for (size_t i = 1; i < nx_ - 1; i++)
	{
		for (size_t j = 1; j < ny_ - 1; j++)
		{
			urhs_[i][j] = dF_dksi(i, j) + dG_deta(i, j);
		}
	}

}

void FlowSolver::computeResidual()
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

void FlowSolver::extrapBoundary()
{
	//1 左边界 do nothing

	//1 上边界
	for (int i = 0; i < nx_; ++i) {
		u_new_[i][ny_ - 1] = u_new_[i][ny_ - 2];
	}
	//1 下边界
	for (int i = 0; i < nx_; ++i) {
		u_new_[i][0] = u_new_[i][1];
	}
	//1 右边界
	for (int j = 0; j < ny_; ++j) {
		u_new_[nx_ - 1][j] = u_new_[nx_ - 2][j];
	}
}

void FlowSolver::writeTecplotFile(size_t timestep) const {
	//if (timestep % GlobalData::ctrlout != 0)
	//	return;
	const auto& x = grid_->getXCoord();
	const auto& y = grid_->getYCoord();
	std::ofstream outFile("../../" + folder_name_ + "/solution_" + std::to_string(timestep) + ".dat");
	if (!outFile.is_open()) {
		throw std::runtime_error("Failed to open file for writing.");
	}
	outFile << "Title=\"Solution Data\"\n";
	//outFile << "Variables=\"X\", \"Y\", \"U\", \"V\", \"U_err\"\n";
	outFile << "Variables=\"X\", \"Y\", \"U\", \"dF/dksi\", \"dG/deta\", \"U_err\"\n";
	outFile << "Zone I = " << nx_ << ", J = " << ny_ << ", F = POINT" << ", STRANDID=1" << ", SOLUTIONTIME=" << timestep << std::endl;
	for (size_t j = 0; j < ny_; ++j) {
		for (size_t i = 0; i < nx_; ++i) {
			double ref = 1.0;
			double err = std::fabs(ref - u_[i][j]);
			//outFile << x[i][j] << " " << y[i][j] << " " << u_[i][j] << " " << v_[i][j] << " " << err << "\n";
			outFile << x[i][j] << " " << y[i][j] << " " << u_[i][j] << " " << dFdksi_[i][j] << " " << dGdeta_[i][j] << " " << err << "\n";
		}
	}
	outFile.close();
	std::cout << "第 " << timestep << " 步结果输出完成;" << std::endl;

}


//Directly discretized----------------------------------------------------------------------------------------------------
double FlowSolver::dau_dx(size_t i, size_t j) {
	double result = 0;
	if (i == 0)
		result = (u_[i + 1][j] - u_[i][j]) / x_step_;
	else if (i == nx_ - 1)
		result = (u_[i][j] - u_[i - 1][j]) / x_step_;
	else
		result = (u_[i + 1][j] - u_[i - 1][j]) / (2 * x_step_);
	return result;
}
double FlowSolver::dbu_dy(size_t i, size_t j) {
	double result = 0;
	if (j == 0)
		result = (u_[i][j + 1] - u_[i][j]) / y_step_;
	else if (j == ny_ - 1)
		result = (u_[i][j] - u_[i][j - 1]) / y_step_;
	else
		result = (u_[i][j + 1] - u_[i][j - 1]) / (2 * y_step_);
	return result;
}
//convective flux-------------------------------------------------------------------------------------------------------
double FlowSolver::F(size_t i, size_t j, const CoefTrans& trans)
{
	double ksix = trans.ksi_x();
	double ksiy = trans.ksi_y();
	double ksit = trans.ksi_t();
	double J = trans.jacob();
	double result = (ksit / J) * u_[i][j] + (ksix / J) * (1 * u_[i][j]) + (ksiy / J) * ((4.0 / 3.0) * u_[i][j]);
	return result;
}
double FlowSolver::G(size_t i, size_t j, const CoefTrans& trans)
{
	double etax = trans.eta_x();
	double etay = trans.eta_y();
	double etat = trans.eta_t();
	double J = trans.jacob();
	double result = (etat / J) * u_[i][j] + (etax / J) * (1 * u_[i][j]) + (etay / J) * ((4.0 / 3.0) * u_[i][j]);
	return result;
}
//Derivative of convective flux------------------------------------------------------------------------------------
double FlowSolver::dF_dksi(size_t i, size_t j)
{
	size_t index = i * ny_ + j;
	const auto& trans = grid_->getCoefTransVec()[index];
	double result = (F(i + 1, j, trans) - F(i - 1, j, trans)) / 2;//这里需要考虑用什么差分算？如果没有通量分裂的话？
	dFdksi_[i][j] = result;
	return result;
}
double FlowSolver::dG_deta(size_t i, size_t j)
{
	size_t index = i * ny_ + j;
	const auto& trans = grid_->getCoefTransVec()[index];
	double result = (G(i, j + 1, trans) - G(i, j - 1, trans)) / 2;
	dGdeta_[i][j] = result;
	return result;
}

//flux splitting------------------------------------------------------------------------------------
double FlowSolver::SW_positive_ksi(double u, double jacobi, double ksi_x, double ksi_y, double ksi_t)
{
	//return 0.5 * (f + fabs(f));
	double f_plus_i = 0.5 * (1 / jacobi * (u * ksi_t + u * ksi_x + (4.0 / 3.0) * u * ksi_y)
		+ std::fabs(1 / jacobi * (u * ksi_t + u * ksi_x + (4.0 / 3.0) * u * ksi_y)));
	return f_plus_i;
}

double FlowSolver::SW_negative_ksi(double u, double jacobi, double ksi_x, double ksi_y, double ksi_t)
{
	//return 0.5 * (f - fabs(f));
	double f_minus_i = 0.5 * (1 / jacobi * (u * ksi_t + u * ksi_x + (4.0 / 3.0) * u * ksi_y)
		- std::abs(1 / jacobi * (u * ksi_t + u * ksi_x + (4.0 / 3.0) * u * ksi_y)));
	return f_minus_i;
}

double FlowSolver::SW_positive_eta(double u, double jacobi, double eta_x, double eta_y, double eta_t)
{
	//return 0.5 * (g + fabs(g));
	double g_plus_j = 0.5 * (1 / jacobi * (u * eta_t + u * eta_x + (4.0 / 3.0) * u * eta_y)
		+ std::fabs(1 / jacobi * (u * eta_t + u * eta_x + (4.0 / 3.0) * u * eta_y)));
	return g_plus_j;
}

double FlowSolver::SW_negative_eta(double u, double jacobi, double eta_x, double eta_y, double eta_t)
{
	//return 0.5 * (g - fabs(g));
	double g_minus_j = 0.5 * (1 / jacobi * (u * eta_t + u * eta_x + (4.0 / 3.0) * u * eta_y)
		- std::fabs(1 / jacobi * (u * eta_t + u * eta_x + (4.0 / 3.0) * u * eta_y)));
	return g_minus_j;
}