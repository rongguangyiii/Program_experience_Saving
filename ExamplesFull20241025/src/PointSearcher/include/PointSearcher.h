
#ifndef POINTSEARCHER_H  
#define POINTSEARCHER_H  
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Fuzzy_sphere.h>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <cmath>
#include <algorithm>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_2 LocPoint;
typedef CGAL::Search_traits_2<K> TreeTraits;
typedef CGAL::Kd_tree<TreeTraits> Kd_tree;
typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;

class PointSearcher {
private:
	std::vector<LocPoint> points;
	std::unordered_map<size_t, size_t> local_to_original;  // 用于局部索引映射原始索引
	std::unique_ptr<Kd_tree> tree;

public:
	PointSearcher(const std::vector<LocPoint>& input_points, const std::vector<size_t>& original_indices);
	void addPoints(const std::vector<LocPoint>& new_pointsVec, const std::vector<size_t>& addIndicesVec);
	void addPoints(const LocPoint& new_point, const size_t addIndex);
	size_t nearest_point_index(const LocPoint& query);
	double nearest_point_distance(const LocPoint& query);
	LocPoint nearest_point_coordinates(const LocPoint& query);
	std::vector<LocPoint> points_within_radius(const LocPoint& query, double radius);
	std::vector<size_t> indices_within_radius(const LocPoint& query, double radius);

private:
	void rebuild_tree();
};
#endif

/*---------------------------------------------------------------------------------------------
//USAGE:

int main() {
	// 定义点集和原始索引，按顺序存储在 vector 中
	std::vector<LocPoint> input_points = {
		LocPoint(0.0, 1.0),  // 索引3  -0
		LocPoint(5.0, 1.0),  // 索引4  -1
		LocPoint(1.0, 2.0),  // 索引13 -2
		LocPoint(2.0, 3.0),  // 索引12 -3
		LocPoint(5.0, 3.0),  // 索引11 -4
		LocPoint(6.0, 4.0),  // 索引5  -5
		LocPoint(1.0, 5.0),  // 索引7  -6
		LocPoint(3.0, 5.0)   // 索引6  -7
	};
	std::vector<size_t> original_indices = { 3, 4, 13, 12, 11,5,7,6 };  // 原始索引

	// 初始化搜索类
	PointSearcher searcher(input_points, original_indices);

	// 添加单个点
	searcher.addPoints(LocPoint(3.0, 3.0), 10);  // 新点
	//添加多个点
	std::vector<LocPoint> addnew_points = {
	LocPoint(0.0, 2.0),  // 索引21  -0
	LocPoint(2.0, 0.0),  // 索引22  -1
	LocPoint(5.0, 5.0)   // 索引23  -1
	};
	std::vector<size_t> addnew_indices = { 21, 22, 23};
	searcher.addPoints(addnew_points, addnew_indices);

	// 查询点
	LocPoint query(1.0, 2.0);

	// 使用1：查找离查询点最近的点
	int nearest_index = searcher.nearest_point_index(query);
	LocPoint nearest_point = searcher.nearest_point_coordinates(query);
	double nearest_distance = searcher.nearest_point_distance(query);

	std::cout << "最近点的索引: " << nearest_index << "\t";
	std::cout << "最近点的坐标: (" << nearest_point.x() << ", " << nearest_point.y() << ")\t";
	std::cout << "距离: " << nearest_distance << std::endl;

	// 使用2：查找在给定半径范围内的所有点
	double radius = 1.0;
	std::vector<LocPoint> points_in_radius = searcher.points_within_radius(query, radius);
	std::vector<size_t> indices_in_radius = searcher.indices_within_radius(query, radius);

	std::cout << "\n在半径 " << radius << " 范围内的点:\n";
	for (int i = 0; i < points_in_radius.size(); ++i) {
		std::cout << "点坐标: (" << points_in_radius[i].x() << ", " << points_in_radius[i].y() << ")\t";
		std::cout << "点的原始索引: " << indices_in_radius[i] << std::endl;
	}

	return 0;
}
---------------------------------------------------------------------------------------------*/