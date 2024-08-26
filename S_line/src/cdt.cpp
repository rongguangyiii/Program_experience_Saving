
#include "cdt.h"

// ����Լ�����ϵĵ�
void TriManager::insertConstraintLine(const std::vector<Coord>& merged_points) {

	// �� Coord ���͵� merged_points ת��Ϊ Point ����
	std::vector<Point> points;
	points.reserve(merged_points.size()); // Ϊ���Ч��Ԥ����ռ�
	for (const auto& coord : merged_points) {
		points.push_back(toCGALPoint(coord));
	}

	// ʹ��ת����ĵ����Լ����
	insertConstraintLine(points);
}

void TriManager::insertConstraintLine(const std::vector<Point>& points) {
	for (size_t i = 0; i < points.size() - 1; ++i) {
		cdt_.insert(points[i], points[i + 1]);
	}
}

// ������Χ����ĵ�
void TriManager::insertBoundaryPoints(const std::vector<Coord>& merged_points) {

	// �� Coord ���͵� merged_points ת��Ϊ Point ����
	std::vector<Point> points;
	points.reserve(merged_points.size()); // Ϊ���Ч��Ԥ����ռ�
	for (const auto& coord : merged_points) {
		points.push_back(toCGALPoint(coord));
	}

	// ʹ��ת����ĵ����Լ����
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

// ��������ʷֵ��ļ������� Tecplot��
void TriManager::outputToTecplot(std::ofstream& out, const std::string& filename) {

	// Tecplot �ļ�ͷ
	out << "TITLE = \""<< filename <<"\"\n";
	out << "VARIABLES = \"X\", \"Y\"\n";
	out << "ZONE T=\"" << filename << "\", N=" << cdt_.number_of_vertices() << ", E=" << cdt_.number_of_faces() << ", F=FEPOINT, ET=TRIANGLE\n";

	// ������ж���
	std::map<CDT::Vertex_handle, int> vertex_index;
	int index = 1;
	for (auto vit = cdt_.finite_vertices_begin(); vit != cdt_.finite_vertices_end(); ++vit) {
		Point p = vit->point();
		out << p.x() << " " << p.y() << "\n";
		vertex_index[vit] = index++;
	}

	// ��������ε����ӹ�ϵ
	for (auto fit = cdt_.finite_faces_begin(); fit != cdt_.finite_faces_end(); ++fit) {
		out << vertex_index[fit->vertex(0)] << " "
			<< vertex_index[fit->vertex(1)] << " "
			<< vertex_index[fit->vertex(2)] << "\n";
	}
	std::cout  << filename << " output complete !\n";

}
