
#include <iostream>
#include "kdTreee.h"

void PointKDTreeManager::removePointsWithinRadius(const std::vector<Coord>& centers, double radius) {
	std::vector<size_t> pointsToRemove;

	for (const auto& center : centers) {
		std::vector<size_t> nearbyPoints = findPointsWithinRadius(center, radius);
		pointsToRemove.insert(pointsToRemove.end(), nearbyPoints.begin(), nearbyPoints.end());
	}

	// 去除重复索引
	std::sort(pointsToRemove.begin(), pointsToRemove.end());
	pointsToRemove.erase(std::unique(pointsToRemove.begin(), pointsToRemove.end()), pointsToRemove.end());

	// 按索引从大到小顺序删除，以防止前面的删除操作影响后面的索引
	std::sort(pointsToRemove.rbegin(), pointsToRemove.rend());
	for (size_t idx : pointsToRemove) {
		pointsList_.erase(pointsList_.begin() + idx);
	}

	// 重新构建kd树，前提是pointsList_中仍然有点
	if (pointsList_.size() > 1) {
		buildKdTree();
	}
	else{
		std::cout << "KdTree size 为 0！\n";
	}
}

void PointKDTreeManager::buildKdTree() {
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	for (const auto& point : pointsList_) {
		if (mode_ == "3D") {
			points->InsertNextPoint(point.x, point.y, point.z);
		}
		else {  // 默认为2D
			points->InsertNextPoint(point.x, point.y, 0.0);
		}
	}

	vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
	polyData->SetPoints(points);

	kDTree_ = vtkSmartPointer<vtkKdTreePointLocator>::New();
	kDTree_->SetDataSet(polyData);
	kDTree_->BuildLocator();
}

std::vector<size_t> PointKDTreeManager::findPointsWithinRadius(const Coord& center, double radius) {
	double Searched_point[3] = { center.x, center.y, center.z };
	vtkSmartPointer<vtkIdList> pointsWithinRadius = vtkSmartPointer<vtkIdList>::New();
	kDTree_->FindPointsWithinRadius(radius, Searched_point, pointsWithinRadius);

	std::vector<size_t> indices;
	for (vtkIdType i = 0; i < pointsWithinRadius->GetNumberOfIds(); ++i) {
		indices.push_back(pointsWithinRadius->GetId(i));
	}

	return indices;
}
