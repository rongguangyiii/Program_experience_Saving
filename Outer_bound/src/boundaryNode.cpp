#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Alpha_shape_2.h>

// 定义CGAL的内核类型和点类型
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
// 定义Delaunay三角化和Alpha形状的相关类型
typedef CGAL::Alpha_shape_vertex_base_2<K>                   Vb;
typedef CGAL::Alpha_shape_face_base_2<K>                     Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>         Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds>               Delaunay;
typedef CGAL::Alpha_shape_2<Delaunay>                        Alpha_shape_2;
typedef Alpha_shape_2::Alpha_shape_edges_iterator            Alpha_shape_edges_iterator;
typedef Alpha_shape_2::Alpha_shape_vertices_iterator         Alpha_shape_vertices_iterator;
typedef Alpha_shape_2::Point_iterator                        Point_iterator;

// 从文件中读取点
std::vector<Point_2> readPointsFromFile(const std::string& filename) 
{
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

// 将点写入文件
void writePointsToFile(const std::vector<Point_2>& points, const std::string& filename) 
{
    std::string outfilename = "D:/ZARAN3/bin/Debug/";
    outfilename = outfilename + "Data_" + filename + ".dat";
    std::ofstream file(outfilename);
    if (!file.is_open()) {
        std::cerr << "错误：无法打开文件 " << filename << " 进行写入。" << std::endl;
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

// 使用Alpha形状算法提取边界点
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

    // 去除重复点
    std::sort(boundaryPoints.begin(), boundaryPoints.end());
    boundaryPoints.erase(std::unique(boundaryPoints.begin(), boundaryPoints.end()), boundaryPoints.end());

    return boundaryPoints;
}

int main() {
    // 定义输入和输出文件名
    std::string inputFilename = "D:/z_test_4code/20240407test/model/points.txt";
    std::string outputFilename = "boundary_points";

    // 从输入文件读取点
    std::vector<Point_2> points = readPointsFromFile(inputFilename);
    if (points.empty()) {
        std::cerr << "错误：输入文件中没有有效的点。" << std::endl;
        return 1;
    }

    // 设置Alpha形状算法的参数
    double alpha = 0.00001; // 根据需要调整该值以获得所需的形状

    // 使用Alpha形状算法提取边界点
    std::vector<Point_2> boundaryPoints = extractBoundaryPoints(points, alpha);

    // 将边界点写入输出文件
    writePointsToFile(boundaryPoints, outputFilename);
    std::cout << "boundaryNode.cpp：边界点已写入至 " << outputFilename << boundaryPoints.size()<< std::endl;

    return 0;
}
