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
	// ����Լ�����ϵĵ�
	void insertConstraintLine(const std::vector<Coord>& points);
	void insertConstraintLine(const std::vector<Point>& points);

	// ������Χ����ĵ�
	void insertBoundaryPoints(const std::vector<Coord>& points);
	void insertBoundaryPoints(const std::vector<Point>& points);

	// ��������ʷֵ��ļ������� Tecplot��
	void outputToTecplot(std::ofstream& out, const std::string& name);

private:
	Point toCGALPoint(const Coord& coord);
	CDT cdt_;
};

/*-----------------------------------------------------------------------
USAGE: �÷�˵��
	   ������Χ���Լ��������cdt;
	   ����Ҫ��ʽ��ʼ�� cdt_����Ϊ�����Զ�ͨ��Ĭ�Ϲ��캯����ʼ����
 ����������������----------------------------------------------------------------
 int main() {
	TriManager tm;

	// ������ֱ�� x = 1 �ϵ�Լ����
	std::vector<Point> constrain_line = {
		Point(1.0, 0.0),
		Point(1.0, 0.5),
		Point(1.0, 1.0),
		Point(1.0, 1.5),
		Point(1.0, 2.0)
	};

	// ����Լ�����ϵĵ�
	tm.insertConstraintLine(constrain_line);

	// ������Χ������ĵ�
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

	// ������Χ��
	tm.insertBoundaryPoints(bound_side);

	// ��������ʷֵ��ļ�
	tm.outputToTecplot("triang20240825.dat");

	return 0;
 }
-----------------------------------------------------------------------*/
