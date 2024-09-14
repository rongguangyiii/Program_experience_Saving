
#include <iostream>
#include "gridGenerate/include/kdTreee.h"

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

Coord PointKDTreeManager::findNearestPoint(const Coord& point) const {
	// 确保KdTree已经构建
	if (!kDTree_) {
		throw std::runtime_error("KdTree is not built yet.");
	}

	double searchPoint[3] = { point.x, point.y, point.z };
	vtkIdType nearestPointId = kDTree_->FindClosestPoint(searchPoint);

	// 根据找到的点ID获取实际坐标
	double nearestPoint[3];
	vtkSmartPointer<vtkPolyData> polyData = vtkPolyData::SafeDownCast(kDTree_->GetDataSet());
	if (!polyData) {
		throw std::runtime_error("DataSet is not of type vtkPolyData.");
	}
	vtkSmartPointer<vtkPoints> points = polyData->GetPoints();
	points->GetPoint(nearestPointId, nearestPoint);

	// 返回最近的点
	return Coord{ nearestPoint[0], nearestPoint[1], nearestPoint[2] };
}

size_t PointKDTreeManager::findNearestPoint(const double x, const double y, const double z) const {
	// 确保KdTree已经构建
	if (!kDTree_) {
		throw std::runtime_error("KdTree is not built yet.");
	}

	double searchPoint[3] = { x,y, z };
	vtkIdType nearestPointId = kDTree_->FindClosestPoint(searchPoint);

	// 返回最近的点索引
	return nearestPointId;
}