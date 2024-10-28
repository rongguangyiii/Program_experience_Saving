#include "PointSearcher/include/PointSearcher.h"

int main() {
	// ����㼯��ԭʼ��������˳��洢�� vector ��
	std::vector<LocPoint> input_points = {
		LocPoint(0.0, 1.0),  // ����3  -0
		LocPoint(5.0, 1.0),  // ����4  -1
		LocPoint(1.0, 2.0),  // ����13 -2
		LocPoint(2.0, 3.0),  // ����12 -3
		LocPoint(5.0, 3.0),  // ����11 -4
		LocPoint(6.0, 4.0),  // ����5  -5
		LocPoint(1.0, 5.0),  // ����7  -6
		LocPoint(3.0, 5.0)   // ����6  -7
	};
	std::vector<size_t> original_indices = { 3, 4, 13, 12, 11,5,7,6 };  // ԭʼ����

	// ��ʼ��������
	PointSearcher searcher(input_points, original_indices);

	// ��ӵ�����
	searcher.addPoints(LocPoint(3.0, 3.0), 10);  // �µ�
	//��Ӷ����
	std::vector<LocPoint> addnew_points = {
	LocPoint(0.0, 2.0),  // ����21  -0
	LocPoint(2.0, 0.0),  // ����22  -1
	LocPoint(5.0, 5.0)   // ����23  -1
	};
	std::vector<size_t> addnew_indices = { 21, 22, 23}; 
	searcher.addPoints(addnew_points, addnew_indices);

	// ��ѯ��
	LocPoint query(1.0, 2.0);

	// ʹ��1���������ѯ������ĵ�
	int nearest_index = searcher.nearest_point_index(query);
	LocPoint nearest_point = searcher.nearest_point_coordinates(query);
	double nearest_distance = searcher.nearest_point_distance(query);

	std::cout << "����������: " << nearest_index << "\t";
	std::cout << "����������: (" << nearest_point.x() << ", " << nearest_point.y() << ")\t";
	std::cout << "����: " << nearest_distance << std::endl;

	// ʹ��2�������ڸ����뾶��Χ�ڵ����е�
	double radius = 1.0;
	std::vector<LocPoint> points_in_radius = searcher.points_within_radius(query, radius);
	std::vector<size_t> indices_in_radius = searcher.indices_within_radius(query, radius);

	std::cout << "\n�ڰ뾶 " << radius << " ��Χ�ڵĵ�:\n";
	for (int i = 0; i < points_in_radius.size(); ++i) {
		std::cout << "������: (" << points_in_radius[i].x() << ", " << points_in_radius[i].y() << ")\t";
		std::cout << "���ԭʼ����: " << indices_in_radius[i] << std::endl;
	}

	return 0;
}
