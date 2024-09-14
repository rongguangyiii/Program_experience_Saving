#pragma once
#include "gridGenerate/include/gridBase.h"
#include "gridGenerate/include/gridBoundS.h"
#include "gridGenerate/include/gridUniform.h"
#include "gridGenerate/include/coefTrans.h"

class CtrlFlow
{
public:
	CtrlFlow():isSpit_(false), isSwallow_(false){};
	~CtrlFlow() {};
	//1.处理边界
	void treatBound();
	std::vector<CoordPtr> foundMoveGridBound();
	//2.标记背景网格点
	void markBackgrid();//可以考虑用kdtree
	void addGrid(std::shared_ptr<GridBase> curgrid) { gridVec_.push_back(std::move(curgrid)); }  // 使用 std::move 转移指针所有权
	void outplot_B(std::ofstream& out, const std::string& zone_title);
	void outplot_M(std::ofstream& out, const std::string& zone_title);
	void solve();
	void move(double vel,double dt);
	void moveGrid(double vel, double dt);
	//bool isthroughput();
	bool isSwallowPoints();
	bool isSpitPoints();
	void calSpitPointsCoesfTrans();
	void foundSpitPoints();
	void genOldGridNeighbors();

	void calTranscoord();
	void calCoefTime(const std::shared_ptr<GridBase>& curgrid);
	void calCoefNoTime(const std::shared_ptr<GridBase>& curgrid);
	//void toflowsolve();
	void markBackWallBound();
	void updateNeighbors();

	//返回gridVec_
	std::vector<std::shared_ptr<GridBase>>& getGridVec() { return gridVec_; }

	//返回allcoefTrans_
	std::unordered_map<PointId, CoefTrans, PointIdHash>& getallcoefTrans() { return allcoefTrans_; }
	std::unordered_map<PointId, CoefTrans, PointIdHash>& getallcoefTrans_old() { return allcoefTrans_old_; }
	std::unordered_map<PointId, CoefTrans, PointIdHash>& getSpitcoefTrans() { return spitoutPointsCoefTrans_; }
	std::vector<CoordPtr>& getSpitoutPoints() { return spitoutPoints_; }
	std::vector<CoordPtr>& getConnectPoints() { return OldconnectPoints_; }
	Coord getoldCoord(PointId index);
	bool isSpit_;
	bool isSwallow_;
private:

	std::vector <std::shared_ptr<GridBase>> gridVec_; //因为 std::unique_ptr 不能被复制，只能通过移动（std::move）来转移所有权。
	//std::shared_ptr<GridBase> moveOldGrid_;
	std::unordered_map<PointId, CoefTrans, PointIdHash> allcoefTrans_;
	std::unordered_map<PointId, CoefTrans, PointIdHash> allcoefTrans_old_;
	std::vector<Coord> movepoints_;
	std::unordered_map<PointId, CoefTrans, PointIdHash> spitoutPointsCoefTrans_;
	std::vector<CoordPtr> spitoutPoints_;
	std::vector<CoordPtr> OldconnectPoints_;
	std::unordered_map<PointId, std::vector<Coord>, PointIdHash> oldmoveGridNeighbors;
	//std::vector<CoordPtr> allPoints_;
};
