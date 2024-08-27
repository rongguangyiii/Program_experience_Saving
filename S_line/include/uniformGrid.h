#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <unordered_map>
#include <set>


struct Coord {
	double x, y, z;
	// z Ĭ��ֵΪ 0.0��ʹ�ýṹ��������ڶ�ά����ά��
	Coord(double x, double y, double z = 0.0) : x(x), y(y), z(z) {}

	bool operator<(const Coord& other) const {
		return std::tie(x, y) < std::tie(other.x, other.y);  // ֻ�Ƚ� x �� y
	}

	bool operator==(const Coord& other) const {
		return std::fabs(x - other.x) < 1e-9 &&
			std::fabs(y - other.y) < 1e-9 &&
			std::fabs(z - other.z) < 1e-9;
	}

};

struct CoordHash {
	std::size_t operator()(const Coord& c) const {
		std::hash<double> double_hash;
		std::size_t h1 = double_hash(c.x);
		std::size_t h2 = double_hash(c.y);
		std::size_t h3 = double_hash(c.z);
		return h1 ^ (h2 << 1) ^ (h3 << 2);
	}
};

class UniformGrid {
public:
	UniformGrid(double x_start, double x_end, double y_start, double y_end, int x_points, int y_points)
		: x_start_(x_start), x_end_(x_end), y_start_(y_start), y_end_(y_end),
		x_points_(x_points), y_points_(y_points) {
		generatePoints();
		generateElements();
	}

	void outputToTecplot(std::ofstream& out, const std::string& zone_title) const;

	std::unordered_map<Coord, std::set<Coord>, CoordHash> getNeighbors() const;

	//����points_
	std::vector<Coord> getPoints() { return points_; }

private:
	double x_start_, x_end_, y_start_, y_end_;
	size_t x_points_, y_points_;
	std::vector<Coord> points_;
	std::vector<std::vector<size_t>> elements_;
	void generatePoints();
	void generateElements();
	void addNeighbor(std::unordered_map<Coord, std::set<Coord>, CoordHash>& map, const Coord& p1, const Coord& p2) const;

};

/*-------------------------------------------------------------------------------
 USAGE: �ṩ���η�Χ����������������
 -------------------------------------------------------------------------------
 int main() {
	// ��һ����������
	UniformGrid grid1(0.0, 4.0, 0.0, 4.0, 21, 21);
	// �ڶ�����������
	UniformGrid grid2(6.0, 10.0, 0.0, 4.0, 21, 21);
	// ������ļ�
	std::ofstream out("uniform_zones_with_elements.dat");
	// �����һ������
	grid1.outputToTecplot("Zone 1", out);
	// ����ڶ�������
	grid2.outputToTecplot("Zone 2", out);
	out.close();
	std::cout << "Points and elements output to uniform_points_zones_with_elements.dat\n";
	return 0;
 }
 -------------------------------------------------------------------------------
 ʾ��2��

int main() {
	// ����һ�� 3x3 ������
	UniformGrid grid(0.0, 2.0, 0.0, 2.0, 3, 3);

	// ��ȡ�ھ���Ϣ
	auto neighbors = grid.getNeighbors();

	// ��ӡÿ������ھ���Ϣ
	for (const auto& [point, neighbor_points] : neighbors) {
		std::cout << "Point (" << point.x << ", " << point.y << ") has neighbors:\n";
		for (const auto& neighbor : neighbor_points) {
			std::cout << "\t(" << neighbor.x << ", " << neighbor.y <<")\n";
		}
	}

	// �������Tecplot�ļ�
	std::ofstream outfile("uniform_grid_output.dat");
	grid.outputToTecplot(outfile, "Uniform Grid");
	outfile.close();

	return 0;
}
--------------------------------------------------------------*/
