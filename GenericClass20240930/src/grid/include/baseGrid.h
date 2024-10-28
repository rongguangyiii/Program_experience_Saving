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
	virtual ~GridBase() = default;  // ������������֤����������ʱ�ܹ���ȷ����
	// ���麯����Ҫ��������ʵ��
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
	//���ز�����=
	virtual std::shared_ptr<GridBase> operator=(const std::shared_ptr<GridBase>& other) = 0;

	// ͨ�÷���
	std::vector<NodePtr>& getPoints() { return pointsVec_; }
	std::vector<Element>& getelements() { return elementsVec_; }
	std::unordered_map<NodeId, std::vector<NodePtr>, NodeIdHash>& getReneighbors() { return neighborsVec_; }
protected:
	std::vector<NodePtr> pointsVec_;
	std::vector<Element> elementsVec_;
	std::unordered_map<NodeId, std::vector<NodePtr>, NodeIdHash> neighborsVec_;
};

/*-----------------------------------------------------------------------------------------------
* USGAE �÷�1
#include <iostream>
#include <fstream>
#include <memory>
#include "gridGenerate/include/BaseGrid.h"
#include "gridGenerate/include/uniformGrid.h"
#include "gridGenerate/include/boundSGrid.h"

int main() {
	// ʹ�ö�̬�Դ�����ͬ���͵��������
	std::unique_ptr<BaseGrid> grid1 = std::make_unique<UniformGrid>(0.0, 1.0, 0.0, 1.0, 11, 11);
	std::unique_ptr<BaseGrid> grid2 = std::make_unique<BoundSGrid>(1.0, 2.0, 0.0, 1.0, 8, 11);

	// ����� Tecplot �ļ�
	std::ofstream fileout("testgrid01.dat");
	if (fileout.is_open()) {
		grid1->outputToTecplot(fileout, "Uniform Grid");
		grid2->outputToTecplot(fileout, "BoundS Grid");
		fileout.close();
	}
	else {
		std::cerr << "Error: Could not open file to write Uniform Grid!" << std::endl;
	}

	//// ʹ�û����ͨ�÷���
	auto& neighbors2 = grid2->getneighbors();

	return 0;
}
-----------------------------------------------------------------------------------------------*/
/*
// * USGAE �÷�2
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
	const std::shared_ptr<GridBase>& cur = Tools::genGrid(0);//uni��0, sin��1, dic��2, perturb��3
	Tools::outGrid(cur, "uni");
	return 0;
}
-----------------------------------------------------------------------------------------------*/
/*
	CoordPtr p1 = std::make_shared<Coord>(1.0, 2.0);
	CoordPtr p2 = std::make_shared<Coord>(3.0, 4.0);
	//��� map ��û�� p1 ��������ֱ�Ӳ���һ���µļ�ֵ�ԡ�����ͨ��ʹ�� map.insert ��ʹ�� map �� operator[] ��������ֵ�ԣ�

	map[p1] = std::set<CoordPtr>();  // ʹ�� operator[] �����¼�ֵ��, ��� p1 �Ѿ����ڣ��������κ����顣
	map.insert({p1, std::set<CoordPtr>()});  // ����ʹ�� insert �����ֵ��
	map[p1].insert(p2); // ����ʹ�� insert �����ֵ��

	map.erase(p1);  // ɾ�� p1 ���������е�ֵ
	map[p1].erase(p2);  // ɾ�� p1 ����Ӧ�ļ����е� p2 ֵ

	if (map.find(p1) != map.end())
	{
	// �� p1 ����
	}

-----------------------------------------------------------------------------------------------*/
