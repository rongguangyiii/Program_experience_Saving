
#include <iostream>
#include "gridGenerate/include/kdTreee.h"

void PointKDTreeManager::removePointsWithinRadius(const std::vector<Coord>& centers, double radius) {
	std::vector<size_t> pointsToRemove;

	for (const auto& center : centers) {
		std::vector<size_t> nearbyPoints = findPointsWithinRadius(center, radius);
		pointsToRemove.insert(pointsToRemove.end(), nearbyPoints.begin(), nearbyPoints.end());
	}

	// ȥ���ظ�����
	std::sort(pointsToRemove.begin(), pointsToRemove.end());
	pointsToRemove.erase(std::unique(pointsToRemove.begin(), pointsToRemove.end()), pointsToRemove.end());

	// �������Ӵ�С˳��ɾ�����Է�ֹǰ���ɾ������Ӱ����������
	std::sort(pointsToRemove.rbegin(), pointsToRemove.rend());
	for (size_t idx : pointsToRemove) {
		pointsList_.erase(pointsList_.begin() + idx);
	}

	// ���¹���kd����ǰ����pointsList_����Ȼ�е�
	if (pointsList_.size() > 1) {
		buildKdTree();
	}
	else{
		std::cout << "KdTree size Ϊ 0��\n";
	}
}

void PointKDTreeManager::buildKdTree() {
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	for (const auto& point : pointsList_) {
		if (mode_ == "3D") {
			points->InsertNextPoint(point.x, point.y, point.z);
		}
		else {  // Ĭ��Ϊ2D
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
	// ȷ��KdTree�Ѿ�����
	if (!kDTree_) {
		throw std::runtime_error("KdTree is not built yet.");
	}

	double searchPoint[3] = { point.x, point.y, point.z };
	vtkIdType nearestPointId = kDTree_->FindClosestPoint(searchPoint);

	// �����ҵ��ĵ�ID��ȡʵ������
	double nearestPoint[3];
	vtkSmartPointer<vtkPolyData> polyData = vtkPolyData::SafeDownCast(kDTree_->GetDataSet());
	if (!polyData) {
		throw std::runtime_error("DataSet is not of type vtkPolyData.");
	}
	vtkSmartPointer<vtkPoints> points = polyData->GetPoints();
	points->GetPoint(nearestPointId, nearestPoint);

	// ��������ĵ�
	return Coord{ nearestPoint[0], nearestPoint[1], nearestPoint[2] };
}

size_t PointKDTreeManager::findNearestPoint(const double x, const double y, const double z) const {
	// ȷ��KdTree�Ѿ�����
	if (!kDTree_) {
		throw std::runtime_error("KdTree is not built yet.");
	}

	double searchPoint[3] = { x,y, z };
	vtkIdType nearestPointId = kDTree_->FindClosestPoint(searchPoint);

	// ��������ĵ�����
	return nearestPointId;
}