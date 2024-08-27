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
	// 插入约束线上的点
	void insertConstraintLine(const std::vector<Coord>& points);
	void insertConstraintLine(const std::vector<Point>& points);

	// 插入外围区域的点
	void insertBoundaryPoints(const std::vector<Coord>& points);
	void insertBoundaryPoints(const std::vector<Point>& points);
	// 返回每个点的所有邻居点
	std::unordered_map<Coord, std::vector<Coord>, CoordHash> getNeighbors();
	// 输出三角剖分到文件（用于 Tecplot）
	void outputToTecplot(std::ofstream& out, const std::string& name);

private:
	Point toCGALPoint(const Coord& coord);
	Coord toCoord(const Point& p) const;
	CDT cdt_;
};

/*-----------------------------------------------------------------------
USAGE: 用法说明
	   给定外围点和约束点生成cdt;
	   不需要显式初始化 cdt_，因为它会自动通过默认构造函数初始化。
 ————————----------------------------------------------------------------
int main() {
    // 创建TriManager对象
    TriManager tm;

    // 定义在直线 x = 1 上的约束点
    std::vector<Coord> constrain_line = {
        {1.0, 0.0, 0.0},
        {1.0, 0.5, 0.0},
        {1.0, 1.0, 0.0},
        {1.0, 1.5, 0.0},
        {1.0, 2.0, 0.0}
    };

    // 插入约束线上的点
    tm.insertConstraintLine(constrain_line);

    // 定义外围区域的点
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

    // 插入外围区域的点
    tm.insertBoundaryPoints(boundary_points);

    // 获取每个点的所有邻居点
    auto neighbors = tm.getNeighbors();

    // 输出邻居点信息
    for (const auto& [point, neighbor_points] : neighbors) {
        std::cout << "Point (" << point.x << ", " << point.y << ", " << point.z << ") has neighbors:\n";
        for (const auto& neighbor : neighbor_points) {
            std::cout << "\t(" << neighbor.x << ", " << neighbor.y << ", " << neighbor.z << ")\n";
        }
    }

    // 输出三角剖分到Tecplot文件
    std::ofstream outfile("triangulation_output.dat");
    tm.outputToTecplot(outfile, "Triangulation Example");
    outfile.close();

    return 0;
}

-----------------------------------------------------------------------*/
