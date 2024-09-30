#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
//#include <CGAL/Delaunay_mesher_2.h>
#include <vector>
#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Constrained_Delaunay_triangulation_2<K> CDT;
typedef CDT::Point Point;


int main() {
	CDT cdt;

	// ������ֱ�� x = 1 �ϵ�Լ����
	std::vector<Point> constrain_line = {
	Point(1.0, 0.0),
	Point(1.0, 0.5),
	Point(1.0, 1.0),
	Point(1.0, 1.5),
	Point(1.0, 2.0)
	};

	for (size_t i = 0; i < constrain_line.size() - 1; ++i) {
		cdt.insert(constrain_line[i], constrain_line[i + 1]);
	}

	// ������Χ������ĵ�
	std::vector<Point> bound_side = {
		Point(0.5, 0.5),
		Point(0.5, 1.5),
		Point(0.2, 1.0),
		Point(0.8, 0.8),
		Point(0.7, 1.3),
		Point(1.5, 0.5),
		Point(1.5, 1.5),
		Point(1.8, 1.0),
		Point(1.2, 0.8),
		Point(1.7, 1.3),
		Point(0.0, 0.0),
		Point(0.0, 0.5),
		Point(0.0, 1.0),
		Point(0.0, 1.5),
		Point(0.0, 2.0),
		Point(0.5, 2.0),
		Point(1.5, 2.0),
		Point(2.0, 2.0),
		Point(2.0, 1.5),
		Point(2.0, 1.0),
		Point(2.0, 0.5),
		Point(2.0, 0.0),
		Point(1.5, 0.0),
		Point(0.5, 0.0),
	};

	// ����Χ����뵽�����ʷ���
	for (const auto& point : bound_side) {
		cdt.insert(point);
	}

	// ��һ���ļ���������� Tecplot
	std::ofstream out("triang20240825.dat");

	// Tecplot �ļ�ͷ
	out << "TITLE = \"Triangulation\"\n";
	out << "VARIABLES = \"X\", \"Y\"\n";
	out << "ZONE T=\"Zone 1\", N=" << cdt.number_of_vertices() << ", E=" << cdt.number_of_faces() << ", F=FEPOINT, ET=TRIANGLE\n";

	// ������ж���
	std::map<CDT::Vertex_handle, int> vertex_index;
	int index = 1;
	for (auto vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit) {
		Point p = vit->point();
		out << p.x() << " " << p.y() << "\n";
		vertex_index[vit] = index++;
	}

	// ��������ε����ӹ�ϵ
	for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
		out << vertex_index[fit->vertex(0)] << " "
			<< vertex_index[fit->vertex(1)] << " "
			<< vertex_index[fit->vertex(2)] << "\n";
	}

	out.close();
	std::cout << "Triangulation output to triangulation.dat\n";

	return 0;
}
