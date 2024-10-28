
#include "PointSearcher/include/PointSearcher.h"
// 内部函数，用于更新KD树
void PointSearcher::rebuild_tree() {
	tree.reset(new Kd_tree(points.begin(), points.end()));
}

PointSearcher::PointSearcher(const std::vector<LocPoint>& input_points, const std::vector<size_t>& original_indices)
	: points(input_points), tree(new Kd_tree(input_points.begin(), input_points.end())) {
	// 构建索引映射
	if (input_points.size() != original_indices.size()) {
		throw std::runtime_error("索引和点数不匹配!");
	}
	for (size_t i = 0; i < original_indices.size(); ++i) {
		local_to_original[i] = original_indices[i];  // 将局部索引与原始索引映射
	}
}

// 添加新的多个点并重建KD树
void PointSearcher::addPoints(const std::vector<LocPoint>& new_pointsVec, const std::vector<size_t>& addIndicesVec) {
	size_t start_index = points.size();  // 记录添加前的点数量
	points.insert(points.end(), new_pointsVec.begin(), new_pointsVec.end());
	for (size_t i = 0; i < new_pointsVec.size(); ++i) {
		local_to_original[start_index + i] = addIndicesVec[i];  // 维护新点的索引映射
	}
	if (points.size() != local_to_original.size()) {
		throw std::runtime_error("索引和点数不匹配!");
	}
	rebuild_tree();  // 重建KD树
}

// 添加单个点并重建KD树
void PointSearcher::addPoints(const LocPoint& new_point, const size_t addIndex) {
	size_t new_index = points.size();  //记录未添加点前的索引
	points.push_back(new_point);
	local_to_original[new_index] = addIndex;  // 维护新点的索引映射
	if (points.size() != local_to_original.size()) {
		throw std::runtime_error("索引和点数不匹配!");
	}
	rebuild_tree();  // 重建KD树
}

// 1. 查找离查询点最近的点的索引
size_t PointSearcher::nearest_point_index(const LocPoint& query) {
	Neighbor_search search(*tree, query, 1);
	Neighbor_search::iterator nearest_it = search.begin();
	LocPoint nearest = nearest_it->first;
	// 找到最近点在 points 中的局部索引
	size_t local_index = std::distance(points.begin(), std::find(points.begin(), points.end(), nearest));
	// 使用 local_to_original 返回原始索引
	auto orig_it = local_to_original.find(local_index);//find用来查找键(key)
	if (orig_it != local_to_original.end()) {
		return orig_it->second;  // 返回原始索引
	}
	return static_cast<size_t>(-1);  // 如果未找到原始索引，返回 -1 或其他错误值
}

// 1. 查找离查询点最近的点的坐标
LocPoint PointSearcher::nearest_point_coordinates(const LocPoint& query) {
	Neighbor_search search(*tree, query, 1);
	return search.begin()->first;
}

// 1. 查找离查询点最近的点的距离
double PointSearcher::nearest_point_distance(const LocPoint& query) {
	Neighbor_search search(*tree, query, 1);
	return std::sqrt(search.begin()->second);
}

// 2. 给定查询半径，返回范围内的点坐标
std::vector<LocPoint> PointSearcher::points_within_radius(const LocPoint& query, double radius) {
	std::vector<LocPoint> results;
	CGAL::Fuzzy_sphere<TreeTraits> fuzzy_sphere(query, radius);
	tree->search(std::back_inserter(results), fuzzy_sphere);
	return results;
}

// 2. 给定查询半径，返回范围内的点索引
std::vector<size_t> PointSearcher::indices_within_radius(const LocPoint& query, double radius) {
	std::vector<LocPoint> results;
	CGAL::Fuzzy_sphere<TreeTraits> fuzzy_sphere(query, radius);
	tree->search(std::back_inserter(results), fuzzy_sphere);

	std::vector<size_t> indices;
	for (const auto& point : results) {
		auto it = std::find(points.begin(), points.end(), point);
		if (it != points.end()) {
			size_t index = std::distance(points.begin(), it);
			// 找到原始索引
			auto orig_it = local_to_original.find(index);//find用来查找键(key)
			if (orig_it != local_to_original.end()) {
				indices.push_back(orig_it->second);  // 使用原始索引
			}
		}
	}
	return indices;
}
