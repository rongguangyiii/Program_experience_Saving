
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
	std::unordered_map<size_t, size_t> local_to_original;  // ���ھֲ�����ӳ��ԭʼ����
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
	// ����㼯��ԭʼ��������˳��洢�� vector ��
	std::vector<LocPoint> input_points = {
		LocPoint(0.0, 1.0),  // ����3  -0
		LocPoint(5.0, 1.0),  // ����4  -1
		LocPoint(1.0, 2.0),  // ����13 -2
		LocPoint(2.0, 3.0),  // ����12 -3
		LocPoint(5.0, 3.0),  // ����11 -4
		LocPoint(6.0, 4.0),  // ����5  -5
		LocPoint(1.0, 5.0),  // ����7  -6
		LocPoint(3.0, 5.0)   // ����6  -7
	};
	std::vector<size_t> original_indices = { 3, 4, 13, 12, 11,5,7,6 };  // ԭʼ����

	// ��ʼ��������
	PointSearcher searcher(input_points, original_indices);

	// ��ӵ�����
	searcher.addPoints(LocPoint(3.0, 3.0), 10);  // �µ�
	//��Ӷ����
	std::vector<LocPoint> addnew_points = {
	LocPoint(0.0, 2.0),  // ����21  -0
	LocPoint(2.0, 0.0),  // ����22  -1
	LocPoint(5.0, 5.0)   // ����23  -1
	};
	std::vector<size_t> addnew_indices = { 21, 22, 23};
	searcher.addPoints(addnew_points, addnew_indices);

	// ��ѯ��
	LocPoint query(1.0, 2.0);

	// ʹ��1���������ѯ������ĵ�
	int nearest_index = searcher.nearest_point_index(query);
	LocPoint nearest_point = searcher.nearest_point_coordinates(query);
	double nearest_distance = searcher.nearest_point_distance(query);

	std::cout << "����������: " << nearest_index << "\t";
	std::cout << "����������: (" << nearest_point.x() << ", " << nearest_point.y() << ")\t";
	std::cout << "����: " << nearest_distance << std::endl;

	// ʹ��2�������ڸ����뾶��Χ�ڵ����е�
	double radius = 1.0;
	std::vector<LocPoint> points_in_radius = searcher.points_within_radius(query, radius);
	std::vector<size_t> indices_in_radius = searcher.indices_within_radius(query, radius);

	std::cout << "\n�ڰ뾶 " << radius << " ��Χ�ڵĵ�:\n";
	for (int i = 0; i < points_in_radius.size(); ++i) {
		std::cout << "������: (" << points_in_radius[i].x() << ", " << points_in_radius[i].y() << ")\t";
		std::cout << "���ԭʼ����: " << indices_in_radius[i] << std::endl;
	}

	return 0;
}
---------------------------------------------------------------------------------------------*/