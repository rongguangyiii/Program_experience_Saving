#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "grid/include/globalData.h"

struct Point {
    double x, y;
};

class EllipsePolygon {
public:
    EllipsePolygon(double a, double b, size_t num_points) : a_(a), b_(b), num_points_(num_points) {
        generatePoints();
    }
    EllipsePolygon(double a, double b, double length) : a_(a), b_(b) {
        num_points_ = (GlobalData::PI * (a_ + b_)) / length;
        generatePoints();
    }

    // �����Tecplot��ʽ��dat�ļ�
    void writeToTecplot(const std::string& filename) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Failed to open file: " << filename << std::endl;
            return;
        }

        // Tecplot header
        file << "TITLE = \"Ellipse\"\n";
        file << "VARIABLES = \"X\", \"Y\"\n";
        file << "ZONE T=\"Ellipse\", N=" << points_.size() << ", E=" << points_.size() << ", F=FEPOINT, ET=LINESEG\n";

        // Output points
        for (const auto& point : points_) {
            file << point.x << " " << point.y << "\n";
        }

        // Output polygon connectivity (simple closed polygon)
        for (size_t i = 0; i < points_.size()-1; ++i) {
            file << (i+ 1) << " " << (i + 2)<< "\n";
        }
		file << points_.size() << " " << "1\n";  // Close the polygon

        file.close();
    }

private:
    double a_, b_;           // ������Ͷ̰���
    size_t num_points_;         // ����εĵ���
    std::vector<Point> points_;

    // ������Բ�ϵĵ�
    void generatePoints() {
        points_.clear();
        double dt = 2.0 * GlobalData::PI / num_points_;
		for (size_t i = 0; i < num_points_; ++i) {
			double t = i * dt;
			double x = 3.0 + a_ * std::cos(t);
			double y = 5.0 + b_ * std::sin(t);
			points_.push_back({ x, y });
        }
    }
};

int main() {
    // ��Բ�ĳ�����Ϊ3���̰���Ϊ2��ʹ��100�������ɶ����
    //EllipsePolygon ellipse(2.0, 1.0, 36);
    EllipsePolygon ellipse(0.2, 0.1, 0.025);

    // �����Tecplot��ʽ���ļ�
    ellipse.writeToTecplot("ellipse.dat");

    std::cout << "Tecplot file generated: ellipse.dat" << std::endl;
    return 0;
}
