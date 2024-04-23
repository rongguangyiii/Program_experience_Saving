#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Alpha_shape_2.h>

// ����CGAL���ں����ͺ͵�����
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
// ����Delaunay���ǻ���Alpha��״���������
typedef CGAL::Alpha_shape_vertex_base_2<K>                   Vb;
typedef CGAL::Alpha_shape_face_base_2<K>                     Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>         Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds>               Delaunay;
typedef CGAL::Alpha_shape_2<Delaunay>                        Alpha_shape_2;
typedef Alpha_shape_2::Alpha_shape_edges_iterator            Alpha_shape_edges_iterator;
typedef Alpha_shape_2::Alpha_shape_vertices_iterator         Alpha_shape_vertices_iterator;
typedef Alpha_shape_2::Point_iterator                        Point_iterator;

// ���ļ��ж�ȡ��
std::vector<Point_2> readPointsFromFile(const std::string& filename) 
{
    std::vector<Point_2> points;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "�����޷����ļ� " << filename << std::endl;
        return points;
    }

    double x, y;
    while (file >> x >> y) {
        points.emplace_back(Point_2(x, y));
    }

    file.close();
    return points;
}

// ����д���ļ�
void writePointsToFile(const std::vector<Point_2>& points, const std::string& filename) 
{
    std::string outfilename = "D:/ZARAN3/bin/Debug/";
    outfilename = outfilename + "Data_" + filename + ".dat";
    std::ofstream file(outfilename);
    if (!file.is_open()) {
        std::cerr << "�����޷����ļ� " << filename << " ����д�롣" << std::endl;
        return;
    }
    file << "TITLE = \"" << filename << "_BoundaryNode" << "\"\n";
    file << "VARIABLES=X, Y \n";

    file << points.size() << std::endl;
    for (const auto& p : points) {
        file << p.x() << " " << p.y() << std::endl;
    }

    file.close();
}

// ʹ��Alpha��״�㷨��ȡ�߽��
std::vector<Point_2> extractBoundaryPoints(const std::vector<Point_2>& points, double alpha) {
    Delaunay dt;
    dt.insert(points.begin(), points.end());

    Alpha_shape_2 as(dt);
    as.set_alpha(alpha);

    std::vector<Point_2> boundaryPoints;
    for (Alpha_shape_edges_iterator it = as.alpha_shape_edges_begin(); it != as.alpha_shape_edges_end(); ++it) {
        Alpha_shape_2::Segment s = as.segment(*it);
        boundaryPoints.push_back(s.vertex(0));
        boundaryPoints.push_back(s.vertex(1));
    }

    // ȥ���ظ���
    std::sort(boundaryPoints.begin(), boundaryPoints.end());
    boundaryPoints.erase(std::unique(boundaryPoints.begin(), boundaryPoints.end()), boundaryPoints.end());

    return boundaryPoints;
}

int main() {
    // �������������ļ���
    std::string inputFilename = "D:/z_test_4code/20240407test/model/points.txt";
    std::string outputFilename = "boundary_points";

    // �������ļ���ȡ��
    std::vector<Point_2> points = readPointsFromFile(inputFilename);
    if (points.empty()) {
        std::cerr << "���������ļ���û����Ч�ĵ㡣" << std::endl;
        return 1;
    }

    // ����Alpha��״�㷨�Ĳ���
    double alpha = 0.00001; // ������Ҫ������ֵ�Ի���������״

    // ʹ��Alpha��״�㷨��ȡ�߽��
    std::vector<Point_2> boundaryPoints = extractBoundaryPoints(points, alpha);

    // ���߽��д������ļ�
    writePointsToFile(boundaryPoints, outputFilename);
    std::cout << "boundaryNode.cpp���߽����д���� " << outputFilename << boundaryPoints.size()<< std::endl;

    return 0;
}
