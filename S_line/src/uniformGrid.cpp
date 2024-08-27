#include "uniformGrid.h"

void UniformGrid::outputToTecplot(std::ofstream& out, const std::string& zone_title) const {
	out << "TITLE = \"" << zone_title << "\"\n";
	out << "VARIABLES = \"X\", \"Y\"\n";
	out << "ZONE T=\"" << zone_title << "\", N=" << points_.size() << ", E=" << elements_.size() << ", F=FEPOINT, ET=QUADRILATERAL\n";

	for (const auto& point : points_) {
		out << std::fixed << std::setprecision(6) << point.x << " " << point.y << "\n";
	}

	for (const auto& elem : elements_) {
		out << elem[0] << " " << elem[1] << " " << elem[2] << " " << elem[3] << "\n";
	}

	std::cout << zone_title << " output complete !\n";
}

void UniformGrid::generatePoints() {
	double x_step = (x_end_ - x_start_) / (x_points_ - 1);
	double y_step = (y_end_ - y_start_) / (y_points_ - 1);

	for (size_t i = 0; i < x_points_; ++i) {
		for (size_t j = 0; j < y_points_; ++j) {
			points_.emplace_back(x_start_ + i * x_step, y_start_ + j * y_step);
		}
	}
}

void UniformGrid::generateElements() {
	for (size_t i = 0; i < x_points_ - 1; ++i) {
		for (size_t j = 0; j < y_points_ - 1; ++j) {
			size_t p1 = i * y_points_ + j + 1;
			size_t p2 = p1 + 1;
			size_t p3 = p1 + y_points_;
			size_t p4 = p3 + 1;
			elements_.push_back({ p1, p3, p4, p2 });
		}
	}
}

// ��ȡÿ������ھӵ�
std::unordered_map<Coord, std::set<Coord>, CoordHash> UniformGrid::getNeighbors() const {
	//ʹ��std::set<Coord>��֤��Ԫ�ص�Ψһ,td::vector<Coord>���ܱ�֤Ԫ�ص�Ψһ
	std::unordered_map<Coord, std::set<Coord>, CoordHash> neighbors_map;

	// ����ÿ����Ԫ�������ھӹ�ϵ
	for (const auto& elem : elements_) {
		Coord p1 = points_.at(elem[0]-1);
		Coord p2 = points_.at(elem[1]-1);
		Coord p3 = points_.at(elem[2]-1);
		Coord p4 = points_.at(elem[3]-1);

		// �����ھӹ�ϵ
		addNeighbor(neighbors_map, p1, p2);
		addNeighbor(neighbors_map, p2, p3);
		addNeighbor(neighbors_map, p3, p4);
		addNeighbor(neighbors_map, p4, p1);

	}

	return neighbors_map;
}

void UniformGrid::addNeighbor(std::unordered_map<Coord, std::set<Coord>, CoordHash>& map, const Coord& p1, const Coord& p2) const {
	if (p1 == p2) return; // ���⽫�Լ���Ϊ�ھ�

	map[p1].insert(p2);
	map[p2].insert(p1);
}