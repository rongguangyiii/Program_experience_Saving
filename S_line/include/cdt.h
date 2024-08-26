#pragma once
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include "uniformGrid.h"


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Constrained_Delaunay_triangulation_2<K> CDT;
typedef CDT::Point Point;

class TriManager {
public:
	// 插入约束线上的点
	void insertConstraintLine(const std::vector<Coord>& points);
	void insertConstraintLine(const std::vector<Point>& points);

	// 插入外围区域的点
	void insertBoundaryPoints(const std::vector<Coord>& points);
	void insertBoundaryPoints(const std::vector<Point>& points);

	// 输出三角剖分到文件（用于 Tecplot）
	void outputToTecplot(std::ofstream& out, const std::string& name);

private:
	Point toCGALPoint(const Coord& coord);
	CDT cdt_;
};

/*-----------------------------------------------------------------------
USAGE: 用法说明
	   给定外围点和约束点生成cdt;
	   不需要显式初始化 cdt_，因为它会自动通过默认构造函数初始化。
 ————————----------------------------------------------------------------
 int main() {
	TriManager tm;

	// 定义在直线 x = 1 上的约束点
	std::vector<Point> constrain_line = {
		Point(1.0, 0.0),
		Point(1.0, 0.5),
		Point(1.0, 1.0),
		Point(1.0, 1.5),
		Point(1.0, 2.0)
	};

	// 插入约束线上的点
	tm.insertConstraintLine(constrain_line);

	// 定义外围侧区域的点
	std::vector<Point> bound_side = {
		Point(0.5, 0.5),
		Point(0.5, 1.5),
		Point(0.2, 1.0),
		Point(0.8, 0.8),
		Point(0.7, 1.3),
		Point(1.5, 0.5),
		Point(1.5, 1.5),
		Point(1.8, 1.0),
		Point(1.2, 0.8),
		Point(1.7, 1.3),
		Point(0.0, 0.0),
		Point(0.0, 0.5),
		Point(0.0, 1.0),
		Point(0.0, 1.5),
		Point(0.0, 2.0),
		Point(0.5, 2.0),
		Point(1.5, 2.0),
		Point(2.0, 2.0),
		Point(2.0, 1.5),
		Point(2.0, 1.0),
		Point(2.0, 0.5),
		Point(2.0, 0.0),
		Point(1.5, 0.0),
		Point(0.5, 0.0),
	};

	// 插入外围点
	tm.insertBoundaryPoints(bound_side);

	// 输出三角剖分到文件
	tm.outputToTecplot("triang20240825.dat");

	return 0;
 }
-----------------------------------------------------------------------*/
