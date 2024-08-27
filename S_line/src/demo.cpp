
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include "uniformGrid.h"
#include "Parabola.h"
#include "cdt.h"
#include "kdTreee.h"

namespace Tools {
	const double EPSILON = 1e-9; // 容差值

	//合并两个vector
	std::vector<Coord> mergeVectors(const std::vector<Coord>& vector1, const std::vector<Coord>& vector2) {
		std::vector<Coord> merged_points;
		merged_points.reserve(vector1.size() + vector2.size());  // 预先分配足够的空间以提高性能
		merged_points.insert(merged_points.end(), vector1.begin(), vector1.end());
		merged_points.insert(merged_points.end(), vector2.begin(), vector2.end());
		return merged_points;
	}

	//将给定数组按y坐标从小到大排序
	void sortPointsByY(std::vector<Coord>& points) {
		std::sort(points.begin(), points.end(), [](const Coord& a, const Coord& b) {
			return a.y < b.y;
			});
	}
};



//int main() {
//	//1.矩形区域
//	UniformGrid grid1(0.0, 4.5, 0.0, 4.0, 46, 41);
//	UniformGrid grid2(5.5, 10.0, 0.0, 4.0, 46, 41);
//	UniformGrid grid_3(4.5, 5.5, 0.0, 4.0, 11, 41);
//	std::vector<Coord> points_1 = grid1.getPoints();
//	std::vector<Coord> points_2 = grid2.getPoints();
//	std::vector<Coord> points_3 = grid_3.getPoints();
//
//	std::ofstream out("newGrid20240826.dat");
//	grid1.outputToTecplot(out, "Zone 1");
//	grid2.outputToTecplot(out, "Zone 2");
//
//	//2.约束分割S线
//	vector<Vector2d> points_down = { Vector2d(5.0, 0.5), Vector2d(4.7, 1.0), Vector2d(5.0, 2.0) };
//	vector<Vector2d> points_up = { Vector2d(5.0, 2.0), Vector2d(5.3, 3.0), Vector2d(5.0, 3.5) };
//	Parabola right_parabola(points_down, 0.1, "right");
//	Parabola left_parabola(points_up, 0.1, "left");
//	std::vector<Coord> constrain_line_down = right_parabola.getPoints();
//	constrain_line_down.pop_back();
//	std::vector<Coord> constrain_line_up = left_parabola.getPoints();
//	std::vector<Coord> constrain_line1, constrain_line2;
//	for (size_t j = 0; j < 5; ++j) {
//		constrain_line1.emplace_back(5.0, 0 + j * 0.1);
//	}
//	for (size_t j = 0; j < 5; ++j) {
//		constrain_line2.emplace_back(5.0, 3.6 + j * 0.1);
//	}
//	//将约束线合并
//	std::vector<Coord> merged_constrain = Tools::mergeVectors(constrain_line1, constrain_line_down);
//	merged_constrain = Tools::mergeVectors(merged_constrain, constrain_line_up);
//	merged_constrain = Tools::mergeVectors(merged_constrain, constrain_line2);
//
//	//3.删除约束线周围的点
//	PointKDTreeManager manager2D(points_3, "2D");
//	manager2D.removePointsWithinRadius(merged_constrain, 0.1);
//	const std::vector<Coord>& remainingPoints2D = manager2D.getPointsList();
//
//	//4.生成中间连接网格
//	TriManager tm;
//	tm.insertConstraintLine(merged_constrain);
//	//4.插入外围点
//	tm.insertBoundaryPoints(remainingPoints2D);
//	tm.outputToTecplot(out, "Zone 3");
//	Tools::outputToTecplot(merged_constrain, out, "Zone 4");//输出线条用于查看
//	std::cout << "Points and elements output to file\n";
//	out.close();
//
//	return 0;
//}

int main() {
	// 创建一个 3x3 的网格
	UniformGrid grid(0.0, 2.0, 0.0, 2.0, 3, 3);

	// 获取邻居信息
	auto neighbors = grid.getNeighbors();

	// 打印每个点的邻居信息
	for (const auto& [point, neighbor_points] : neighbors) {
		std::cout << "Point (" << point.x << ", " << point.y << ") has neighbors:\n";
		for (const auto& neighbor : neighbor_points) {
			std::cout << "\t(" << neighbor.x << ", " << neighbor.y <<")\n";
		}
	}

	// 输出网格到Tecplot文件
	std::ofstream outfile("uniform_grid_output.dat");
	grid.outputToTecplot(outfile, "Uniform Grid");
	outfile.close();

	return 0;
}