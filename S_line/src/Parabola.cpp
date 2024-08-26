
#include "Parabola.h"

using namespace Eigen;
using namespace std;


// 构造函数：左开口抛物线
Parabola::Parabola(const vector<Vector2d>& points, double step, std::string which) :step_(step), tag_(which) {
	if (tag_ == "left") {
		calLeftParabolaCoefficients(points);
	}
	else if (tag_ == "right") {
		calRightParabolaCoefficients(points);
	}
	else
	{
		std::cout << "the direction is wrong!\n";
	}
	generatePoints();
}

// 生成抛物线上的点
void Parabola::generatePoints()
{
	double start = points_[0][1];
	double end = points_[2][1];
	double epsilon = 5E-6;
	std::function<double(double)> computeX;

	// 根据tag_设置计算x的函数
	if (tag_ == "left") {
		computeX = [this](double y) { return -a_ * y * y + b_ * y + c_; };
	}
	else if (tag_ == "right") {
		computeX = [this](double y) { return a_ * y * y + b_ * y + c_; };
	}

	// 生成数据点
	for (double y = start; y < end + epsilon; y += step_) {
		double x = computeX(y);
		data_points_.emplace_back(x, y);
	}
}


// 输出抛物线到文件
void Parabola::outputToTecplot(ofstream& out, const string& name) const 
{
	out << "TITLE = \"" << name << "\"\n";
	out << "VARIABLES = \"X\", \"Y\"\n";
	out << "ZONE T=\"" << name << "\", N=" << data_points_.size() << ", E=" << data_points_.size() - 1 << ", F=FEPOINT, ET=LINESEG\n";

	for (const auto& point : data_points_) {
		out << point.x << " " << point.y << "\n";
	}

	for (size_t i = 0; i < data_points_.size() - 1; ++i) {
		out << i + 1 << " " << i + 2 << "\n";
	}
	std::cout << name << " output complete !\n";
}

// 计算左开口抛物线的系数
void Parabola::calLeftParabolaCoefficients(const vector<Vector2d>& points) {
	points_ = points;
	MatrixXd A(3, 3);
	VectorXd b(3);

	for (size_t i = 0; i < points.size(); ++i) {
		double y = points[i][1];
		A(i, 0) = -y * y;
		A(i, 1) = y;
		A(i, 2) = 1.0;
		b(i) = points[i][0];
	}

	VectorXd coefficients = A.colPivHouseholderQr().solve(b);
	a_ = coefficients(0);
	b_ = coefficients(1);
	c_ = coefficients(2);
}

// 计算右开口抛物线的系数
void Parabola::calRightParabolaCoefficients(const vector<Vector2d>& points) {
	points_ = points;
	MatrixXd A(3, 3);
	VectorXd b(3);

	for (size_t i = 0; i < points.size(); ++i) {
		double y = points[i][1];
		A(i, 0) = y * y;
		A(i, 1) = y;
		A(i, 2) = 1.0;
		b(i) = points[i][0];
	}

	VectorXd coefficients = A.colPivHouseholderQr().solve(b);
	a_ = coefficients(0);
	b_ = coefficients(1);
	c_ = coefficients(2);
}


