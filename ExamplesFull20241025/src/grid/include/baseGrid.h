#pragma once
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <fstream>
#include <iostream>
#include "grid/include/Coord.h"
#include "grid/include/node.h"
#include "grid/include/element.h"
#include "grid/include/globalData.h"
#include "grid/include/coefTrans.h"

class GridBase {
public:
	GridBase() = default;
	virtual ~GridBase() = default;  // 虚析构函数保证派生类析构时能够正确调用
	// 纯虚函数，要求派生类实现
	virtual void generatePoints() = 0;
	virtual void generateElements() = 0;
	virtual void coefTrans() = 0;
	virtual void outputToTecplot(std::ofstream& out, const std::string& zone_title) const = 0;
	virtual size_t getXnum() const = 0;
	virtual size_t getYnum() const = 0;
	virtual double getXstart() const = 0;
	virtual double getXend() const = 0;
	virtual double getYstart() const = 0;
	virtual double getYend() const = 0;
	virtual std::vector<std::vector<double>> getXCoord() const = 0;
	virtual std::vector<std::vector<double>> getYCoord() const = 0;
	virtual void genNeiborNode() = 0;
	virtual const std::vector<CoefTrans>& getCoefTransVec() const = 0;
	virtual std::shared_ptr<GridBase> operator=(const std::shared_ptr<GridBase>& other) = 0;

	// 通用方法
	std::vector<NodePtr>& getPoints() { return pointsVec_; }
	std::vector<Element>& getelements() { return elementsVec_; }
	std::unordered_map<NodeId, std::vector<NodePtr>, NodeIdHash>& getReneighbors() { return neighborsVec_; }
protected:
	std::vector<NodePtr> pointsVec_;
	std::vector<Element> elementsVec_;
	std::unordered_map<NodeId, std::vector<NodePtr>, NodeIdHash> neighborsVec_;
};

/*-----------------------------------------------------------------------------------------------
* USGAE 用法1
#include <iostream>
#include <fstream>
#include <memory>
#include "gridGenerate/include/BaseGrid.h"
#include "gridGenerate/include/uniformGrid.h"
#include "gridGenerate/include/boundSGrid.h"

int main() {
	// 使用多态性创建不同类型的网格对象
	std::unique_ptr<BaseGrid> grid1 = std::make_unique<UniformGrid>(0.0, 1.0, 0.0, 1.0, 11, 11);
	std::unique_ptr<BaseGrid> grid2 = std::make_unique<BoundSGrid>(1.0, 2.0, 0.0, 1.0, 8, 11);

	// 输出到 Tecplot 文件
	std::ofstream fileout("testgrid01.dat");
	if (fileout.is_open()) {
		grid1->outputToTecplot(fileout, "Uniform Grid");
		grid2->outputToTecplot(fileout, "BoundS Grid");
		fileout.close();
	}
	else {
		std::cerr << "Error: Could not open file to write Uniform Grid!" << std::endl;
	}

	//// 使用基类的通用方法
	auto& neighbors2 = grid2->getneighbors();

	return 0;
}
-----------------------------------------------------------------------------------------------*/
/*
// * USGAE 用法2
#include <iostream>
#include <fstream>
#include <memory>
#include "grid/include/Basegrid.h"
#include "grid/include/gridSin.h"
#include "grid/include/gridDisc.h"
#include "grid/include/gridPerturb.h"
#include "grid/include/gridUniform.h"
#include "grid/include/globalData.h"
#include "Tools/include/tools.h"

int main() {
	const std::shared_ptr<GridBase>& cur = Tools::genGrid(0);//uni：0, sin：1, dic：2, perturb：3
	Tools::outGrid(cur, "uni");
	return 0;
}
-----------------------------------------------------------------------------------------------*/
/*
	CoordPtr p1 = std::make_shared<Coord>(1.0, 2.0);
	CoordPtr p2 = std::make_shared<Coord>(3.0, 4.0);
	//如果 map 中没有 p1 键，可以直接插入一个新的键值对。可以通过使用 map.insert 或使用 map 的 operator[] 来新增键值对：

	map[p1] = std::set<CoordPtr>();  // 使用 operator[] 插入新键值对, 如果 p1 已经存在，不会做任何事情。
	map.insert({p1, std::set<CoordPtr>()});  // 或者使用 insert 插入键值对
	map[p1].insert(p2); // 或者使用 insert 插入键值对

	map.erase(p1);  // 删除 p1 键及其所有的值
	map[p1].erase(p2);  // 删除 p1 键对应的集合中的 p2 值

	if (map.find(p1) != map.end())
	{
	// 键 p1 存在
	}

-----------------------------------------------------------------------------------------------*/
