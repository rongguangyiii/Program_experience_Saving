#pragma once
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <fstream>
#include <iostream>
#include "gridGenerate/include/coord.h"
#include "flowSolve/include/globalData.h"
#include "gridGenerate/include/coefTrans.h"


class GridBase {
public:
    GridBase() = default;
    virtual ~GridBase() = default;  // ������������֤����������ʱ�ܹ���ȷ����

    // ���麯����Ҫ��������ʵ��
    virtual void generatePoints() = 0;
    virtual void generateElements() = 0;
    virtual void coordTrans() = 0;
    virtual void outputToTecplot(std::ofstream& out, const std::string& zone_title) const = 0;
    virtual size_t getXnum() const = 0;
    virtual size_t getYnum() const = 0;
    virtual double getXstart() const = 0;
    virtual double getXend() const = 0;
    virtual double getYstart() const = 0;
    virtual double getYend() const = 0;
    virtual void genNeiborNode() = 0;
   //���ز�����=
    virtual std::shared_ptr<GridBase> operator=(const std::shared_ptr<GridBase>& other) = 0;

	// ͨ�÷���
    void setCoefTrans(const CoordPtr& p, const CoefTrans& coef);
    std::vector<CoordPtr>& getPoints() { return points_; }
    //CoordPtr& getPoints(PointId index);
    std::vector<Element>& getelements() { return elements_; }
    std::unordered_map<PointId, std::vector<CoordPtr>, PointIdHash>& getReneighbors() { return Reneighbors_; }
    std::unordered_map<PointId, CoefTrans, PointIdHash>& getTransCoef() { return coefTrans_; }
protected:
    std::vector<CoordPtr> points_;
    std::vector<Element> elements_;
    std::unordered_map<PointId, std::vector<CoordPtr>, PointIdHash> Reneighbors_;
    std::unordered_map<PointId, CoefTrans, PointIdHash> coefTrans_;
};

/*-----------------------------------------------------------------------------------------------
* USGAE �÷�ʾ��
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

*/