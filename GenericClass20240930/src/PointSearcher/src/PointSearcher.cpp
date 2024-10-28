
#include "PointSearcher/include/PointSearcher.h"
// �ڲ����������ڸ���KD��
void PointSearcher::rebuild_tree() {
	tree.reset(new Kd_tree(points.begin(), points.end()));
}

PointSearcher::PointSearcher(const std::vector<LocPoint>& input_points, const std::vector<size_t>& original_indices)
	: points(input_points), tree(new Kd_tree(input_points.begin(), input_points.end())) {
	// ��������ӳ��
	if (input_points.size() != original_indices.size()) {
		throw std::runtime_error("�����͵�����ƥ��!");
	}
	for (size_t i = 0; i < original_indices.size(); ++i) {
		local_to_original[i] = original_indices[i];  // ���ֲ�������ԭʼ����ӳ��
	}
}

// ����µĶ���㲢�ؽ�KD��
void PointSearcher::addPoints(const std::vector<LocPoint>& new_pointsVec, const std::vector<size_t>& addIndicesVec) {
	size_t start_index = points.size();  // ��¼���ǰ�ĵ�����
	points.insert(points.end(), new_pointsVec.begin(), new_pointsVec.end());
	for (size_t i = 0; i < new_pointsVec.size(); ++i) {
		local_to_original[start_index + i] = addIndicesVec[i];  // ά���µ������ӳ��
	}
	if (points.size() != local_to_original.size()) {
		throw std::runtime_error("�����͵�����ƥ��!");
	}
	rebuild_tree();  // �ؽ�KD��
}

// ��ӵ����㲢�ؽ�KD��
void PointSearcher::addPoints(const LocPoint& new_point, const size_t addIndex) {
	size_t new_index = points.size();  //��¼δ��ӵ�ǰ������
	points.push_back(new_point);
	local_to_original[new_index] = addIndex;  // ά���µ������ӳ��
	if (points.size() != local_to_original.size()) {
		throw std::runtime_error("�����͵�����ƥ��!");
	}
	rebuild_tree();  // �ؽ�KD��
}

// 1. �������ѯ������ĵ������
size_t PointSearcher::nearest_point_index(const LocPoint& query) {
	Neighbor_search search(*tree, query, 1);
	Neighbor_search::iterator nearest_it = search.begin();
	LocPoint nearest = nearest_it->first;
	// �ҵ�������� points �еľֲ�����
	size_t local_index = std::distance(points.begin(), std::find(points.begin(), points.end(), nearest));
	// ʹ�� local_to_original ����ԭʼ����
	auto orig_it = local_to_original.find(local_index);//find�������Ҽ�(key)
	if (orig_it != local_to_original.end()) {
		return orig_it->second;  // ����ԭʼ����
	}
	return static_cast<size_t>(-1);  // ���δ�ҵ�ԭʼ���������� -1 ����������ֵ
}

// 1. �������ѯ������ĵ������
LocPoint PointSearcher::nearest_point_coordinates(const LocPoint& query) {
	Neighbor_search search(*tree, query, 1);
	return search.begin()->first;
}

// 1. �������ѯ������ĵ�ľ���
double PointSearcher::nearest_point_distance(const LocPoint& query) {
	Neighbor_search search(*tree, query, 1);
	return std::sqrt(search.begin()->second);
}

// 2. ������ѯ�뾶�����ط�Χ�ڵĵ�����
std::vector<LocPoint> PointSearcher::points_within_radius(const LocPoint& query, double radius) {
	std::vector<LocPoint> results;
	CGAL::Fuzzy_sphere<TreeTraits> fuzzy_sphere(query, radius);
	tree->search(std::back_inserter(results), fuzzy_sphere);
	return results;
}

// 2. ������ѯ�뾶�����ط�Χ�ڵĵ�����
std::vector<size_t> PointSearcher::indices_within_radius(const LocPoint& query, double radius) {
	std::vector<LocPoint> results;
	CGAL::Fuzzy_sphere<TreeTraits> fuzzy_sphere(query, radius);
	tree->search(std::back_inserter(results), fuzzy_sphere);

	std::vector<size_t> indices;
	for (const auto& point : results) {
		auto it = std::find(points.begin(), points.end(), point);
		if (it != points.end()) {
			size_t index = std::distance(points.begin(), it);
			// �ҵ�ԭʼ����
			auto orig_it = local_to_original.find(index);//find�������Ҽ�(key)
			if (orig_it != local_to_original.end()) {
				indices.push_back(orig_it->second);  // ʹ��ԭʼ����
			}
		}
	}
	return indices;
}
