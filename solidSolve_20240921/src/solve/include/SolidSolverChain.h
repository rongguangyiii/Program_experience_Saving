#pragma once
#include <vector>
#include <memory>
#include "grid/include/baseGrid.h"

class SolidSolverChain
{
public:
	SolidSolverChain(std::shared_ptr<GridBase> g, const size_t nx, const size_t ny);
	void solve();
	void initialize();
	void initialize_ref();
	void computeRHS();
	void computeRHSTrans();
	void solveEuler();
	void solveEulerTrans();
	void extrapBoundary();
	void timeAdvance();
	void computeResidual();
	void calPartialDer();
	void calPartialDer_ref(); 
	void writeTecplotFile(size_t timestep)const;
private:
	double Q(size_t i, size_t j);//du_dx
	double R(size_t i, size_t j);//du_dy
	double S(size_t i, size_t j);//dv_dx
	double T(size_t i, size_t j);//dv_dy
	double du_dksi(size_t i, size_t j);
	double du_deta(size_t i, size_t j);
	double dv_dksi(size_t i, size_t j);
	double dv_deta(size_t i, size_t j);
	double Fu(size_t i, size_t j, const CoefTrans& trans);
	double Fv(size_t i, size_t j, const CoefTrans& trans);
	double Gu(size_t i, size_t j, const CoefTrans& trans);
	double Gv(size_t i, size_t j, const CoefTrans& trans);
	double dFu_dksi(size_t i, size_t j);
	double dFv_deta(size_t i, size_t j);
	double dGu_dksi(size_t i, size_t j);
	double dGv_deta(size_t i, size_t j);
private:
	double FuHalfPlus(size_t i, size_t j, const CoefTrans& trans);
	double FuHalfMius(size_t i, size_t j, const CoefTrans& trans);
	double FvHalfPlus(size_t i, size_t j, const CoefTrans& trans);
	double FvHalfMius(size_t i, size_t j, const CoefTrans& trans);
	double GuHalfPlus(size_t i, size_t j, const CoefTrans& trans);
	double GuHalfMius(size_t i, size_t j, const CoefTrans& trans);
	double GvHalfPlus(size_t i, size_t j, const CoefTrans& trans);
	double GvHalfMius(size_t i, size_t j, const CoefTrans& trans);
private:
	bool useRef_;
	double caltime_, residual_;
	double E_, miu_, dt_;
	double x_step_, y_step_;
	size_t nx_, ny_;
	double c1_, c2_, c3_;
	std::shared_ptr<GridBase> grid_;
	std::vector<std::vector<double>> u_;
	std::vector<std::vector<double>> u_new_;
	std::vector<std::vector<double>> v_new_;
	std::vector<std::vector<double>> v_;
	std::vector<std::vector<double>> urhs_;
	std::vector<std::vector<double>> vrhs_;

	//debug
	std::vector<std::vector<double>> Q_;
	std::vector<std::vector<double>> R_;
	std::vector<std::vector<double>> S_;
	std::vector<std::vector<double>> T_;
	std::vector<std::vector<double>> DFuDksi_;
	std::vector<std::vector<double>> DFvDeta_;
};

