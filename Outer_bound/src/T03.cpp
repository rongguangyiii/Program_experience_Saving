#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Alpha_shape_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef K::Segment_2 Segment_2;
typedef K::FT FT;
typedef CGAL::Alpha_shape_vertex_base_2<K>                   Vb;
typedef CGAL::Alpha_shape_face_base_2<K>                     Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>         Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds>               Delaunay;
typedef CGAL::Alpha_shape_2<Delaunay>                        Alpha_shape_2;
typedef Alpha_shape_2::Alpha_shape_edges_iterator            Alpha_shape_edges_iterator;
typedef Alpha_shape_2::Alpha_shape_vertices_iterator         Alpha_shape_vertices_iterator;
typedef Alpha_shape_2::Point_iterator                        Point_iterator;

std::vector<Point_2> readPointsFromFile(const std::string& filename) {
    std::vector<Point_2> points;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "错误：无法打开文件 " << filename << std::endl;
        return points;
    }

    double x, y;
    while (file >> x >> y) {
        points.emplace_back(Point_2(x, y));
    }

    file.close();
    return points;
}

void writePointsToFile(const std::vector<Point_2>& points, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "错误：无法打开文件 " << filename << " 进行写入。" << std::endl;
        return;
    }

    for (const auto& p : points) {
        file << p.x() << " " << p.y() << std::endl;
    }

    file.close();
}

std::vector<Point_2> extractBoundaryPoints(const std::vector<Point_2>& points, double alpha) {
    Delaunay dt;
    dt.insert(points.begin(), points.end());

    Alpha_shape_2 as(dt);
    as.set_alpha(alpha);

    std::vector<Segment_2> segments;
    for (Alpha_shape_edges_iterator it = as.alpha_shape_edges_begin(); it != as.alpha_shape_edges_end(); ++it) {
        segments.push_back(as.segment(*it));
    }

    std::vector<Point_2> boundaryPoints;
    for (const auto& segment : segments) {
        boundaryPoints.push_back(segment.source());
        boundaryPoints.push_back(segment.target());
    }

    // Remove duplicate points
    std::sort(boundaryPoints.begin(), boundaryPoints.end());
    boundaryPoints.erase(std::unique(boundaryPoints.begin(), boundaryPoints.end()), boundaryPoints.end());

    return boundaryPoints;
}

int main() {
    std::string inputFilename = "D:/z_test_4code/20240407test/model/points.txt";
    std::string outputFilename = "boundary_points.txt";

    std::vector<Point_2> points = readPointsFromFile(inputFilename);
    if (points.empty()) {
        std::cerr << "错误：输入文件中没有有效的点。" << std::endl;
        return 1;
    }

    double alpha = 0.00001;

    std::vector<Point_2> boundaryPoints = extractBoundaryPoints(points, alpha);

    writePointsToFile(boundaryPoints, outputFilename);
    std::cout << "有"<< boundaryPoints.size() <<"个边界点已写入至 " << outputFilename << std::endl;

    return 0;
}
