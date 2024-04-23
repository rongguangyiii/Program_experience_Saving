#include <iostream>  
#include <vector>  
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>  
#include <CGAL/Delaunay_triangulation_2.h>  
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>  
#include <CGAL/Constrained_Delaunay_triangulation_2.h>  
#include <CGAL/Polygon_2.h>  
#include <CGAL/Polygon_with_holes_2.h>  

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef CGAL::Polygon_with_holes_2<K> Polygon_with_holes_2;

int main() {
    // ��������poly��hole�Ķ���  
    std::vector<Point_2> poly_vertices = { Point_2(1, 1), Point_2(3, 1), Point_2(5, 1), Point_2(5, 3), Point_2(6, 5), Point_2(3, 6), Point_2(1, 5), Point_2(2, 3) };
    std::vector<Point_2> hole_vertices = { Point_2(3, 2), Point_2(4, 2), Point_2(4, 4), Point_2(3, 4) };

    // ��������κʹ����Ķ����  
    Polygon_2 poly(poly_vertices.begin(), poly_vertices.end());
    Polygon_with_holes_2 poly_with_hole(poly, hole_vertices.begin(), hole_vertices.end());

    // ����Լ��Delaunay�����ʷֶ���  
    CGAL::Constrained_Delaunay_triangulation_2<K> cdt;

    // �������εı���ΪԼ��  
    cdt.insert_constraint(poly.edges_begin(), poly.edges_end());

    // ����hole�ı���ΪԼ��  
    cdt.insert_constraint(hole_vertices.begin(), hole_vertices.end());

    // ��������ʷֵĽ��  
    for (CGAL::Constrained_Delaunay_triangulation_2<K>::Finite_faces_iterator it = cdt.finite_faces_begin(); it != cdt.finite_faces_end(); ++it) {
        CGAL::Triangle_2<K> triangle = cdt.triangle(*it);
        std::cout << "Triangle: " << triangle[0] << ", " << triangle[1] << ", " << triangle[2] << std::endl;
    }

    return 0;
}