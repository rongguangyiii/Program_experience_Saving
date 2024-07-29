#include"Coordinate/include/Coordinate.h"
#include <Eigen/Dense>
//#include <iostream>
#include <vector>
#include <cmath>
#include <map>

// 定义一个结构体来存储点的坐标和各种属性
struct PointData {
	double temperature_;
	double density_;
	double pressure_;
};

//反距离插值
class IDWInterpolator {
public:
	// 构造函数，初始化已知点的坐标和温度
	IDWInterpolator(const std::vector<Coordinate>& coords, const std::vector<PointData>& temps, double p)
		: coordinates_(coords), pointData_(temps), power_(p) {}

	// 计算IDW插值
	std::map<std::string, double> interpolate(const Coordinate& targetP) const;

private:
	std::vector<Coordinate> coordinates_;
	std::vector<PointData> pointData_;
	double power_;

};

// 克里金插值协方差模型(半变异函数模型)
enum class VariogramModel {
	LINEAR,
	SPHERICAL,
	EXPONENTIAL,
	GAUSSIAN
};
// 克里金插值
class KrigingInterpolator {
public:
	KrigingInterpolator(const std::vector<Coordinate>& points, const std::vector<PointData>& data, VariogramModel model = VariogramModel::LINEAR, double a = 0, double b = 0)
		: points_(points), data_(data), model_(model), cofa_(a), cofb_(b) {}

	// 计算Kriging插值
	std::map<std::string, double> interpolate(const Coordinate& target);

private:
	std::vector<std::pair<double, double>> experimentalVariogram(const std::vector<Coordinate>& points, const std::vector<double>& values);
	std::pair<double, double> fitVariogramModel(const std::vector<std::pair<double, double>>& variogram);
	double variogramFunction(double h)const;
	double interpolateValue(const std::vector<double>& values, const Coordinate& target);  // 添加这一行

	std::vector<Coordinate> points_;
	std::vector<PointData> data_;
	VariogramModel model_;
	double cofa_;
	double cofb_;
};
