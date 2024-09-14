#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <unordered_map>
#include <set>
#include "gridGenerate/include/coord.h"

namespace GlobalData
{
	const std::string meshType = "Uniform";//网格类型 Uniform; Semicircular;
	const std::string calmethod = "Common";//计算RHS方法 Deer(冻结系数); Common(使用计算点度量系数); Basic(无通量分裂)
	const size_t nx = 101;
	const size_t ny = 51;
	const double cfl = 0.9;      // CFL数
	const double ctrltime = 100;      // 总时间
	const size_t ctrlstep = 10000;      // 时间步数
	const double tolerance = 1e-16;  // 收敛阈值
	const size_t ctrlout = 50; //输出控制
	const size_t movestep = 1; 
	const double PI = 3.141592653589793238462643383279502884197169399;//圆周率pi

	void createFolder(const std::string& basePath);
	bool isPointInPolygon_moreInner(const CoordPtr& p, const std::vector<CoordPtr>& polygon);

    // 计算质心
	Coord calculateCentroid(const std::vector<CoordPtr>& points);
	//Coord calculateCentroid(const std::set<CoordPtr>& points);

    // 计算给定点相对于质心的角度（0 到 2π）
	double calAngle_02pi(const Coord& centroid, const CoordPtr& point);
	double calAngle_pipi(const Coord& centroid, const CoordPtr& point);

    // 按逆时针方向（从0到2π）排序
	void sortPointsCounterclockwise(std::vector<CoordPtr>& points);

	double angleDifference(double angle1, double angle2);

	std::vector<CoordPtr> filterPointsByAngleRange(const std::vector<CoordPtr>& points, double angleStart, double angleEnd);

	CoordPtr& PointAtAngle(std::vector<CoordPtr>& points, double targetAngle);
}

namespace Tools
{
	bool isCollinear(const Coord& p1, const Coord& p2, const Coord& p3);
	bool isPointInPolygon_moreOuter(const Coord& p, const std::vector<Coord>& polygon);
	bool isPointInPolygon_moreInner(const Coord& p, const std::vector<Coord>& polygon);
	void reorderPointsIndex(std::vector<std::pair<Coord, size_t>>& PointLists, const Coord& modelCenter);
	double Distance(const CoordPtr& point1, const CoordPtr& point2);
	std::vector<CoordPtr> cood2CoodPtr(const std::vector<Coord>& points);
	std::vector<Coord> coodPtr2Cood(const std::vector<CoordPtr>& pointsPtr);
	// 通用检查函数
	template<typename T, typename K>
	const T& checkKeyExists(const std::unordered_map<K, T, PointIdHash>& map, const K& key, const std::string& mapName) {
		auto it = map.find(key);
		if (it == map.end()){
			std::cout << "Key not found in" << mapName << "\n";
			throw std::runtime_error("Key not found!!!");
		}
		return it->second;
	}
}