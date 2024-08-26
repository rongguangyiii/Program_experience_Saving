
#include <iostream>
#include "kdTreee.h"

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
