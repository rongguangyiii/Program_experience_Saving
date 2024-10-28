#pragma once
#include <vector>
#include <memory>
#include "grid/include/baseGrid.h"
#include "solve/include/gradient.h"

class FlowSolveND2UP
{
public:
	FlowSolveND2UP(std::shared_ptr<GridBase> g, const size_t nx, const size_t ny);
	void solve();
	void initialize();
	void initialize_ref();
	void computeRHSTrans();
	void solveEulerTrans();
	void extrapBoundary();
	void timeAdvance();
	void computeResidual();
private:
	void writeTecplotFile(size_t timestep)const;
	double F(size_t i, size_t j);
	double G(size_t i, size_t j);
	double dF_dksi(size_t i, size_t j);
	double dG_deta(size_t i, size_t j);
	double getRef(const double x, const double y) const;
	std::pair<double, double> calculateErrorNorms() const;
	// flux splitting
	double SW_positive_ksi(double u, double jacobi, double ksi_x, double ksi_y, double ksi_t);
	double SW_negative_ksi(double u, double jacobi, double ksi_x, double ksi_y, double ksi_t);
	double SW_positive_eta(double u, double jacobi, double eta_x, double eta_y, double eta_t);
	double SW_negative_eta(double u, double jacobi, double eta_x, double eta_y, double eta_t);
private:
	bool useRef_;
	double caltime_, residual_;
	double x_step_, y_step_, dt_;
	size_t nx_, ny_;
	std::string folder_name_;
	std::shared_ptr<GridBase> grid_;
	std::vector<std::vector<double>> u_;
	std::vector<std::vector<double>> u_new_;
	std::vector<std::vector<double>> urhs_;

	//debug
	std::vector<std::vector<double>> Q_;
	std::vector<std::vector<double>> R_;
	std::vector<std::vector<double>> S_;
	std::vector<std::vector<double>> T_;
	std::vector<std::vector<double>> dFdksi_;
	std::vector<std::vector<double>> dGdeta_;
};
