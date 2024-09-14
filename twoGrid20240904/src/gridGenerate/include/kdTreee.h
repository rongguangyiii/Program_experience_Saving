#pragma once
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkKdTree.h>
#include <vtkKdTreePointLocator.h>
#include <vtkIdList.h>
#include <vector>
#include <string>
#include <algorithm>
#include "gridGenerate/include/coord.h"

class PointKDTreeManager {
public:
	PointKDTreeManager(const std::vector<Coord>& pointsList, const std::string& mode = "2D")
		: pointsList_(pointsList), mode_(mode) {
		if (mode_ != "2D" && mode_ != "3D") {
			throw std::invalid_argument("Invalid mode: must be '2D' or '3D'");
		}
		if (!pointsList_.empty()) {
			buildKdTree();
		}
	}

	void removePointsWithinRadius(const std::vector<Coord>& centers, double radius);
	std::vector<size_t> findPointsWithinRadius(const Coord& center, double radius);
	const std::vector<Coord>& getPointsList() const { return pointsList_; }
	Coord findNearestPoint(const Coord& point) const;
	size_t findNearestPoint(const double x, const double y, const double z) const;
private:
	std::vector<Coord> pointsList_;
	vtkSmartPointer<vtkKdTreePointLocator> kDTree_;
	std::string mode_;  // "2D" 或 "3D"

	void buildKdTree();
	//std::vector<size_t> findPointsWithinRadius(const Coord& center, double radius);
};


/*---------------------------------------------------------------------------------

 USAGE: 用法示例
		给定点数据，和控制点数据，删除控制点给定半径内的点
 ---------------------------------------------------------------------------------
 int main() {
	// 定义一些二维点
	std::vector<Coord> points2D = { {0, 0}, {1, 1}, {2, 2}, {3, 3}, {4, 4}, {5, 5} , {1, 4} };
	PointKDTreeManager manager2D(points2D, "2D");

	std::vector<Coord> centers2D = { {1, 1}, {3, 3} };
	double radius2D = 1.5;

	manager2D.removePointsWithinRadius(centers2D, radius2D);

	const std::vector<Coord>& remainingPoints2D = manager2D.getPointsList();
	std::cout << "Remaining 2D Points: " << std::endl;
	for (const auto& point : remainingPoints2D) {
		std::cout << "(" << point.x << ", " << point.y << ")" << std::endl;
	}

	// 定义一些三维点
	std::vector<Coord> points3D = { {0, 0, 0}, {1, 1, 1}, {2, 2, 2}, {3, 3, 3}, {4, 4, 4} };
	PointKDTreeManager manager3D(points3D, "3D");

	std::vector<Coord> centers3D = { {1, 1, 1}, {3, 3, 3} };
	double radius3D = 1.8;

	manager3D.removePointsWithinRadius(centers3D, radius3D);

	const std::vector<Coord>& remainingPoints3D = manager3D.getPointsList();
	std::cout << "Remaining 3D Points: " << std::endl;
	for (const auto& point : remainingPoints3D) {
		std::cout << "(" << point.x << ", " << point.y << ", " << point.z << ")" << std::endl;
	}

	return 0;
 }

---------------------------------------------------------------------------------*/