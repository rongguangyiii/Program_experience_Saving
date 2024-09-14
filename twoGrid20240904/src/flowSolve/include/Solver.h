#pragma once
#include <vector>
#include <string>
#include <memory>
//#include"grid/include/grid.h"
#include "gridGenerate/include/coord.h"
#include "flowSolve/include/globalData.h"
//#include "gridGenerate/include/coefTrans.h"
#include "flowSolve/include/ctrlFlow.h"
#include "gridGenerate/include/gridBase.h"

class FlowSolver
{
public:
	FlowSolver(std::shared_ptr<CtrlFlow> &ctrlflow) : ctrlflow_(ctrlflow), dt(0.0025),
													  caltime(0.0), residual(1.0), isfirst_(true){}
	FlowSolver() = delete;
	~FlowSolver() {};

	void solve();
	void outplot_B(std::ofstream& out, const std::string& zone_title, const size_t timestep);
	void outplot_M(std::ofstream& out, const std::string& zone_title, const size_t timestep);
private:
	void initialflow();
	void FilterPointsToCal();
	void FilterPointsToCal(const bool);
	void saveOldNeighbors();

	void timeAdvance();
	void calSpitPoints();
	void flowBoundary();
	void extrapBoundary_b();
	void extrapBoundary_s();
	void solverEuler();
	void computeResidual();
	void computeSpitPointsRHS_deer();
	void computeRHS_deer(); //冻结系数
	void computeRHS_common();//各自度量系数
	void writeTecplotFile(size_t timestep);
	// 通量分裂函数
	double SW_positive_ksi(double u, double jacobi, double ksi_x, double ksi_y, double ksi_t);
	double SW_negative_ksi(double u, double jacobi, double ksi_x, double ksi_y, double ksi_t);
	double SW_positive_eta(double u, double jacobi, double eta_x, double eta_y, double eta_t);
	double SW_negative_eta(double u, double jacobi, double eta_x, double eta_y, double eta_t);

	double caltime;
	double dt;
	double residual;
	bool isfirst_;
	std::shared_ptr<CtrlFlow>& ctrlflow_;
	std::unordered_map<PointId, double, PointIdHash> u;
	std::unordered_map<PointId, double, PointIdHash> u_new;
	std::unordered_map<PointId, double, PointIdHash> rhs;
	std::unordered_map<PointId, double, PointIdHash> spit_rhs;
	std::vector<CoordPtr> calPoints_;
	std::unordered_map<PointId, std::vector<CoordPtr>, PointIdHash> AllNeighbors_;
	std::unordered_map<PointId, std::vector<Coord>, PointIdHash> OldNeighbors_;


};



