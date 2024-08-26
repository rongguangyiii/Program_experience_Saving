
#include "cdt.h"

// 插入约束线上的点
void TriManager::insertConstraintLine(const std::vector<Coord>& merged_points) {

	// 将 Coord 类型的 merged_points 转换为 Point 类型
	std::vector<Point> points;
	points.reserve(merged_points.size()); // 为提高效率预分配空间
	for (const auto& coord : merged_points) {
		points.push_back(toCGALPoint(coord));
	}

	// 使用转换后的点插入约束线
	insertConstraintLine(points);
}

void TriManager::insertConstraintLine(const std::vector<Point>& points) {
	for (size_t i = 0; i < points.size() - 1; ++i) {
		cdt_.insert(points[i], points[i + 1]);
	}
}

// 插入外围区域的点
void TriManager::insertBoundaryPoints(const std::vector<Coord>& merged_points) {

	// 将 Coord 类型的 merged_points 转换为 Point 类型
	std::vector<Point> points;
	points.reserve(merged_points.size()); // 为提高效率预分配空间
	for (const auto& coord : merged_points) {
		points.push_back(toCGALPoint(coord));
	}

	// 使用转换后的点插入约束线
	insertBoundaryPoints(points);
}

void TriManager::insertBoundaryPoints(const std::vector<Point>& points) {
	for (const auto& point : points) {
		cdt_.insert(point);
	}
}

Point TriManager::toCGALPoint(const Coord& coord) {
	return Point(coord.x, coord.y);
}

// 输出三角剖分到文件（用于 Tecplot）
void TriManager::outputToTecplot(std::ofstream& out, const std::string& filename) {

	// Tecplot 文件头
	out << "TITLE = \""<< filename <<"\"\n";
	out << "VARIABLES = \"X\", \"Y\"\n";
	out << "ZONE T=\"" << filename << "\", N=" << cdt_.number_of_vertices() << ", E=" << cdt_.number_of_faces() << ", F=FEPOINT, ET=TRIANGLE\n";

	// 输出所有顶点
	std::map<CDT::Vertex_handle, int> vertex_index;
	int index = 1;
	for (auto vit = cdt_.finite_vertices_begin(); vit != cdt_.finite_vertices_end(); ++vit) {
		Point p = vit->point();
		out << p.x() << " " << p.y() << "\n";
		vertex_index[vit] = index++;
	}

	// 输出三角形的连接关系
	for (auto fit = cdt_.finite_faces_begin(); fit != cdt_.finite_faces_end(); ++fit) {
		out << vertex_index[fit->vertex(0)] << " "
			<< vertex_index[fit->vertex(1)] << " "
			<< vertex_index[fit->vertex(2)] << "\n";
	}
	std::cout  << filename << " output complete !\n";

}
