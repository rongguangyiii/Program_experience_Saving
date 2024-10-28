#include "grid/include/globalData.h"
#include "solve/include/flowSolveND1UP.h"
#include "tools/include/tools.h"
#include <cmath>

FlowSolveND1UP::FlowSolveND1UP(const std::shared_ptr<GridBase> g, size_t nx, size_t ny) : residual_(0.0), grid_(g), caltime_(0),
u_(nx, std::vector<double>(ny, 0.0)), u_new_(nx, std::vector<double>(ny, 0.0)),
urhs_(nx, std::vector<double>(ny, 0.0)), Q_(nx, std::vector<double>(ny, 0.0)), R_(nx, std::vector<double>(ny, 0.0)),
S_(nx, std::vector<double>(ny, 0.0)), T_(nx, std::vector<double>(ny, 0.0)),
dFdksi_(nx, std::vector<double>(ny, 0.0)), dGdeta_(nx, std::vector<double>(ny, 0.0))
{
	folder_name_ = "resultNoDEER_1ordUP";
	double x_start_ = grid_->getXstart();
	double y_start_ = grid_->getYstart();
	double x_end_ = grid_->getXend();
	double y_end_ = grid_->getYend();
	nx_ = nx;
	ny_ = ny;
	x_step_ = (x_end_ - x_start_) / (nx - 1);
	y_step_ = (y_end_ - y_start_) / (ny - 1);
	double c1 =4/3.0;
	double c = 0.15;//控制系数，控制时间步长小于稳定性条件。
	dt_ = c * x_step_ / c1;
	useRef_ = true;
};

void FlowSolveND1UP::solve()
{
	Tools::createFolder("../../"+ folder_name_);
	if (useRef_) {
		initialize_ref();//使用理论值初始化场内变量,边界不用外推
	}
	else {
		initialize();
	}
	timeAdvance();
}

void FlowSolveND1UP::initialize() {
	const auto& Ycoord = grid_->getYCoord();
	const auto& Xcoord = grid_->getXCoord();
	for (size_t i = 0; i < nx_; i++)
	{
		for (size_t j = 0; j < ny_; j++)
		{
			u_[i][j] = 0.0;
		}
	}
	//左边界 u=0固定;
	for (int j = 0; j < ny_; ++j) {
		//u_[0][j] = 1.0;
		u_[0][j] = 1.0 + Ycoord[0][j] - (4.0 / 3.0) * Xcoord[0][j];
	}
	//下边界
	for (int i = 0; i < nx_; ++i) {
		u_[i][0] = 1.0 + Ycoord[i][0] - (4.0 / 3.0) * Xcoord[i][0];
	}

	u_new_ = u_;
}

void FlowSolveND1UP::initialize_ref() {
	const auto& Xcoord = grid_->getXCoord();
	const auto& Ycoord = grid_->getYCoord();
	for (size_t i = 0; i < nx_; i++)
	{
		for (size_t j = 0; j < ny_; j++)
		{
			u_[i][j] = 1.0 + Ycoord[i][j] - (4.0 / 3.0) * Xcoord[i][j];
		}
	}
	u_new_ = u_;
}

void FlowSolveND1UP::timeAdvance()
{
	size_t timestep = 0;
	residual_ = 1.0;
	writeTecplotFile(timestep);
	while (residual_ > GlobalData::tolerance)
	//while (timestep < 50)
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

void FlowSolveND1UP::solveEulerTrans()
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
	if (!useRef_)
	{
		extrapBoundary();
	}
	computeResidual();
}

void FlowSolveND1UP::computeRHSTrans()
{
	for (size_t i = 1; i < nx_ - 1; i++)
	{
		for (size_t j = 1; j < ny_ - 1; j++)
		{
			urhs_[i][j] = dF_dksi(i, j) + dG_deta(i, j);
		}
	}

}

void FlowSolveND1UP::computeResidual()
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

void FlowSolveND1UP::extrapBoundary()
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

void FlowSolveND1UP::writeTecplotFile(size_t timestep) const {
	if (timestep % GlobalData::ctrlout != 0)
		return;
	const auto& x = grid_->getXCoord();
	const auto& y = grid_->getYCoord();
	std::ofstream outFile("../../" + folder_name_ + "/solution_" + std::to_string(timestep) + ".dat");
	if (!outFile.is_open()) {
		std::cerr << "Failed to open file for writing." << std::endl;
		throw std::runtime_error("Failed to open file for writing.");
	}
	outFile << "Title=\"Solution Data\"\n";
	//outFile << "Variables=\"X\", \"Y\", \"U\", \"V\", \"U_err\"\n";
	outFile << "Variables=\"X\", \"Y\", \"U\", \"dF/dksi\", \"dG/deta\", \"U_err\"\n";
	outFile << "Zone I = " << nx_ << ", J = " << ny_ << ", F = POINT" << ", STRANDID=1" << ", SOLUTIONTIME=" << timestep << std::endl;
	double err_l2_norm = 0.0;
	double err_inf_norm = 0.0;
	double ref = 10e9;

	for (size_t j = 0; j < ny_; ++j) {
		for (size_t i = 0; i < nx_; ++i) {
			ref = getRef(x[i][j], y[i][j]);
			double err = std::fabs(ref - u_[i][j]);
			outFile << x[i][j] << " " << y[i][j] << " " << u_[i][j] << " " << dFdksi_[i][j] << " " << dGdeta_[i][j] << " " << err << "\n";
		}
	}

	std::tie(err_l2_norm, err_inf_norm) = calculateErrorNorms();
	outFile << "第 " << timestep << " 步结果输出完成; L2范数: " << err_l2_norm << ", 无穷范数: " << err_inf_norm << std::endl;

	outFile.close();
}
double FlowSolveND1UP::getRef(const double x, const double y)const
{
	//double ref = 1.0;
	double ref = 1.0 + y - (4.0 / 3.0) * x;
	return ref;
}

std::pair<double, double> FlowSolveND1UP::calculateErrorNorms() const {
	double err_l2_norm = 0.0;
	double err_inf_norm = 0.0;
	double ref = 10e9;
	const auto& x = grid_->getXCoord();
	const auto& y = grid_->getYCoord();

	for (size_t j = 0; j < ny_; ++j) {
		for (size_t i = 0; i < nx_; ++i) {
			ref = getRef(x[i][j], y[i][j]);
			double err = std::fabs(ref - u_[i][j]);
			err_l2_norm += err * err;
			if (err > err_inf_norm) {
				err_inf_norm = err;
			}
		}
	}

	err_l2_norm = std::sqrt(err_l2_norm);
	err_l2_norm /= (nx_ * ny_);
	return { err_l2_norm, err_inf_norm };
}

//convective flux-------------------------------------------------------------------------------------------------------
double FlowSolveND1UP::F(size_t i, size_t j)
{
	size_t index = i * ny_ + j;
	const CoefTrans& trans = grid_->getCoefTransVec()[index];
	double ksix = trans.ksi_x();
	double ksiy = trans.ksi_y();
	double ksit = trans.ksi_t();
	double J = trans.jacob();
	double result = (ksit / J) * u_[i][j] + (ksix / J) * (1 * u_[i][j]) + (ksiy / J) * ((4.0 / 3.0) * u_[i][j]);
	return result;
}
double FlowSolveND1UP::G(size_t i, size_t j)
{
	size_t index = i * ny_ + j;
	const CoefTrans& trans = grid_->getCoefTransVec()[index];
	double etax = trans.eta_x();
	double etay = trans.eta_y();
	double etat = trans.eta_t();
	double J = trans.jacob();
	double result = (etat / J) * u_[i][j] + (etax / J) * (1 * u_[i][j]) + (etay / J) * ((4.0 / 3.0) * u_[i][j]);
	return result;
}
//Derivative of convective flux------------------------------------------------------------------------------------
double FlowSolveND1UP::dF_dksi(size_t i, size_t j)
{
	//double result = (F(i + 1, j) - F(i - 1, j)) / 2;//二阶中心: f(i+1)-f(i-1)/2ΔH
	double result = (F(i, j) - F(i - 1, j));//一阶迎风前差:f(i)+f(i-1)/ΔH，后差:f(i+1)-f(i)/ΔH
	dFdksi_[i][j] = result;
	return result;
}
double FlowSolveND1UP::dG_deta(size_t i, size_t j)
{
	//double result = (G(i, j + 1) - G(i, j - 1)) / 2;
	double result = (G(i, j) - G(i, j - 1));
	dGdeta_[i][j] = result;
	return result;
}

//flux splitting------------------------------------------------------------------------------------
double FlowSolveND1UP::SW_positive_ksi(double u, double jacobi, double ksi_x, double ksi_y, double ksi_t)
{
	//return 0.5 * (f + fabs(f));
	double f_plus_i = 0.5 * (1 / jacobi * (u * ksi_t + u * ksi_x + (4.0 / 3.0) * u * ksi_y))
		+ std::fabs(1 / jacobi * (u * ksi_t + u * ksi_x + (4.0 / 3.0) * u * ksi_y));
	return f_plus_i;
}
double FlowSolveND1UP::SW_negative_ksi(double u, double jacobi, double ksi_x, double ksi_y, double ksi_t)
{
	//return 0.5 * (f - fabs(f));
	double f_minus_i = 0.5 * (1 / jacobi * (u * ksi_t + u * ksi_x + (4.0 / 3.0) * u * ksi_y))
		- std::abs(1 / jacobi * (u * ksi_t + u * ksi_x + (4.0 / 3.0) * u * ksi_y));
	return f_minus_i;
}
double FlowSolveND1UP::SW_positive_eta(double u, double jacobi, double eta_x, double eta_y, double eta_t)
{
	//return 0.5 * (g + fabs(g));
	double g_plus_j = 0.5 * (1 / jacobi * (u * eta_t + u * eta_x + (4.0 / 3.0) * u * eta_y))
		+ std::fabs(1 / jacobi * (u * eta_t + u * eta_x + (4.0 / 3.0) * u * eta_y));
	return g_plus_j;
}
double FlowSolveND1UP::SW_negative_eta(double u, double jacobi, double eta_x, double eta_y, double eta_t)
{
	//return 0.5 * (g - fabs(g));
	double g_minus_j = 0.5 * (1 / jacobi * (u * eta_t + u * eta_x + (4.0 / 3.0) * u * eta_y))
		- std::fabs(1 / jacobi * (u * eta_t + u * eta_x + (4.0 / 3.0) * u * eta_y));
	return g_minus_j;
}
