#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>

struct Coord {
    double x, y,z;
    // z 默认值为 0.0，使得结构体可以用于二维和三维点
	Coord(double x, double y, double z = 0.0) : x(x), y(y), z(z) {}
};

class UniformGrid {
public:
    UniformGrid(double x_start, double x_end, double y_start, double y_end, int x_points, int y_points)
        : x_start_(x_start), x_end_(x_end), y_start_(y_start), y_end_(y_end),
        x_points_(x_points), y_points_(y_points) {
        generatePoints();
        generateElements();
    }

    void outputToTecplot(std::ofstream& out, const std::string& zone_title) const;
   
    //返回points_
    std::vector<Coord> getPoints() { return points_; }

private:
    double x_start_, x_end_, y_start_, y_end_;
    size_t x_points_, y_points_;
    std::vector<Coord> points_;
    std::vector<std::vector<size_t>> elements_;
    void generatePoints();

    void generateElements();
};

/*-------------------------------------------------------------------------------
 USAGE: 提供矩形范围的输入来生成网格
 -------------------------------------------------------------------------------
 int main() {
    // 第一个矩形区域
    UniformGrid grid1(0.0, 4.0, 0.0, 4.0, 21, 21);

    // 第二个矩形区域
    UniformGrid grid2(6.0, 10.0, 0.0, 4.0, 21, 21);

    // 打开输出文件
    std::ofstream out("uniform_zones_with_elements.dat");

    // 输出文件头
    out << "TITLE = \"Uniform Points with Elements\"\n";
    out << "VARIABLES = \"X\", \"Y\"\n";

    // 输出第一个区域
    grid1.outputToTecplot("Zone 1", out);

    // 输出第二个区域
    grid2.outputToTecplot("Zone 2", out);

    out.close();
    std::cout << "Points and elements output to uniform_points_zones_with_elements.dat\n";

    return 0;
 } 
--------------------------------------------------------------*/
