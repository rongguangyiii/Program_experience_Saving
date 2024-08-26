
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include "uniformGrid.h"
#include "Parabola.h"
#include "cdt.h"
#include "kdTreee.h"

namespace Tools {
	//�ϲ�����vector
	std::vector<Coord> mergeVectors(const std::vector<Coord>& vector1, const std::vector<Coord>& vector2) {
		std::vector<Coord> merged_points;
		merged_points.reserve(vector1.size() + vector2.size());  // Ԥ�ȷ����㹻�Ŀռ����������
		merged_points.insert(merged_points.end(), vector1.begin(), vector1.end());
		merged_points.insert(merged_points.end(), vector2.begin(), vector2.end());
		return merged_points;
	}

	//���������鰴y�����С��������
	void sortPointsByY(std::vector<Coord>& points) {
		std::sort(points.begin(), points.end(), [](const Coord& a, const Coord& b) {
			return a.y < b.y;
			});
	}


	void outputToTecplot(const std::vector<Coord>& pointslist, ofstream& out, const string& name)
	{
		out << "TITLE = \"" << name << "\"\n";
		out << "VARIABLES = \"X\", \"Y\"\n";
		out << "ZONE T=\"" << name << "\", N=" << pointslist.size() << ", E=" << pointslist.size() - 1 << ", F=FEPOINT, ET=LINESEG\n";

		for (const auto& point : pointslist) {
			out << point.x << " " << point.y << "\n";
		}

		for (size_t i = 0; i < pointslist.size() - 1; ++i) {
			out << i + 1 << " " << i + 2 << "\n";
		}

		std::cout << name << " output complete !\n";
	}
}

int main() {
	//1.��������
	UniformGrid grid1(0.0, 4.5, 0.0, 4.0, 46, 41);
	UniformGrid grid2(5.5, 10.0, 0.0, 4.0, 46, 41);
	UniformGrid grid_3(4.5, 5.5, 0.0, 4.0, 11, 41);
	std::vector<Coord> points_1 = grid1.getPoints();
	std::vector<Coord> points_2 = grid2.getPoints();
	std::vector<Coord> points_3 = grid_3.getPoints();

	std::ofstream out("newGrid20240826.dat");
	grid1.outputToTecplot(out, "Zone 1");
	grid2.outputToTecplot(out, "Zone 2");

	//2.Լ���ָ�S��
	vector<Vector2d> points_down = { Vector2d(5.0, 0.5), Vector2d(4.7, 1.0), Vector2d(5.0, 2.0) };
	vector<Vector2d> points_up = { Vector2d(5.0, 2.0), Vector2d(5.3, 3.0), Vector2d(5.0, 3.5) };
	Parabola right_parabola(points_down, 0.1, "right");
	Parabola left_parabola(points_up, 0.1, "left");
	std::vector<Coord> constrain_line_down = right_parabola.getPoints();
	constrain_line_down.pop_back();
	std::vector<Coord> constrain_line_up = left_parabola.getPoints();
	std::vector<Coord> constrain_line1, constrain_line2;
	for (size_t j = 0; j < 5; ++j) {
		constrain_line1.emplace_back(5.0, 0 + j * 0.1);
	}
	for (size_t j = 0; j < 5; ++j) {
		constrain_line2.emplace_back(5.0, 3.6 + j * 0.1);
	}
	//��Լ���ߺϲ�
	std::vector<Coord> merged_constrain = Tools::mergeVectors(constrain_line1, constrain_line_down);
	merged_constrain = Tools::mergeVectors(merged_constrain, constrain_line_up);
	merged_constrain = Tools::mergeVectors(merged_constrain, constrain_line2);

	//3.ɾ��Լ������Χ�ĵ�
	PointKDTreeManager manager2D(points_3, "2D");
	manager2D.removePointsWithinRadius(merged_constrain, 0.1);
	const std::vector<Coord>& remainingPoints2D = manager2D.getPointsList();

	//4.�����м���������
	TriManager tm;
	tm.insertConstraintLine(merged_constrain);
	//4.������Χ��
	tm.insertBoundaryPoints(remainingPoints2D);
	tm.outputToTecplot(out, "Zone 3");
	Tools::outputToTecplot(merged_constrain, out, "Zone 4");//����������ڲ鿴
	std::cout << "Points and elements output to file\n";
	out.close();

	return 0;
}
