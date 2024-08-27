#pragma once
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include "uniformGrid.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> TDS;
typedef CGAL::Exact_predicates_tag Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag> CDT;
typedef CDT::Point Point;
typedef CDT::Vertex_handle Vertex_handle;

class TriManager {
public:
	// ����Լ�����ϵĵ�
	void insertConstraintLine(const std::vector<Coord>& points);
	void insertConstraintLine(const std::vector<Point>& points);

	// ������Χ����ĵ�
	void insertBoundaryPoints(const std::vector<Coord>& points);
	void insertBoundaryPoints(const std::vector<Point>& points);
	// ����ÿ����������ھӵ�
	std::unordered_map<Coord, std::vector<Coord>, CoordHash> getNeighbors();
	// ��������ʷֵ��ļ������� Tecplot��
	void outputToTecplot(std::ofstream& out, const std::string& name);

private:
	Point toCGALPoint(const Coord& coord);
	Coord toCoord(const Point& p) const;
	CDT cdt_;
};

/*-----------------------------------------------------------------------
USAGE: �÷�˵��
	   ������Χ���Լ��������cdt;
	   ����Ҫ��ʽ��ʼ�� cdt_����Ϊ�����Զ�ͨ��Ĭ�Ϲ��캯����ʼ����
 ����������������----------------------------------------------------------------
int main() {
    // ����TriManager����
    TriManager tm;

    // ������ֱ�� x = 1 �ϵ�Լ����
    std::vector<Coord> constrain_line = {
        {1.0, 0.0, 0.0},
        {1.0, 0.5, 0.0},
        {1.0, 1.0, 0.0},
        {1.0, 1.5, 0.0},
        {1.0, 2.0, 0.0}
    };

    // ����Լ�����ϵĵ�
    tm.insertConstraintLine(constrain_line);

    // ������Χ����ĵ�
    std::vector<Coord> boundary_points = {
        {0.5, 0.5, 0.0},
        {0.5, 1.5, 0.0},
        {0.2, 1.0, 0.0},
        {0.8, 0.8, 0.0},
        {0.7, 1.3, 0.0},
        {1.5, 0.5, 0.0},
        {1.5, 1.5, 0.0},
        {1.8, 1.0, 0.0},
        {1.2, 0.8, 0.0},
        {1.7, 1.3, 0.0},
        {0.0, 0.0, 0.0},
        {0.0, 0.5, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 1.5, 0.0},
        {0.0, 2.0, 0.0},
        {0.5, 2.0, 0.0},
        {1.5, 2.0, 0.0},
        {2.0, 2.0, 0.0},
        {2.0, 1.5, 0.0},
        {2.0, 1.0, 0.0},
        {2.0, 0.5, 0.0},
        {2.0, 0.0, 0.0},
        {1.5, 0.0, 0.0},
        {0.5, 0.0, 0.0},
    };

    // ������Χ����ĵ�
    tm.insertBoundaryPoints(boundary_points);

    // ��ȡÿ����������ھӵ�
    auto neighbors = tm.getNeighbors();

    // ����ھӵ���Ϣ
    for (const auto& [point, neighbor_points] : neighbors) {
        std::cout << "Point (" << point.x << ", " << point.y << ", " << point.z << ") has neighbors:\n";
        for (const auto& neighbor : neighbor_points) {
            std::cout << "\t(" << neighbor.x << ", " << neighbor.y << ", " << neighbor.z << ")\n";
        }
    }

    // ��������ʷֵ�Tecplot�ļ�
    std::ofstream outfile("triangulation_output.dat");
    tm.outputToTecplot(outfile, "Triangulation Example");
    outfile.close();

    return 0;
}

-----------------------------------------------------------------------*/
