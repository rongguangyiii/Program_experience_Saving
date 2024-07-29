#include"Coordinate/include/Coordinate.h"
#include <Eigen/Dense>
//#include <iostream>
#include <vector>
#include <cmath>
#include <map>

// ����һ���ṹ�����洢�������͸�������
struct PointData {
	double temperature_;
	double density_;
	double pressure_;
};

//�������ֵ
class IDWInterpolator {
public:
	// ���캯������ʼ����֪���������¶�
	IDWInterpolator(const std::vector<Coordinate>& coords, const std::vector<PointData>& temps, double p)
		: coordinates_(coords), pointData_(temps), power_(p) {}

	// ����IDW��ֵ
	std::map<std::string, double> interpolate(const Coordinate& targetP) const;

private:
	std::vector<Coordinate> coordinates_;
	std::vector<PointData> pointData_;
	double power_;

};

// ������ֵЭ����ģ��(����캯��ģ��)
enum class VariogramModel {
	LINEAR,
	SPHERICAL,
	EXPONENTIAL,
	GAUSSIAN
};
// ������ֵ
class KrigingInterpolator {
public:
	KrigingInterpolator(const std::vector<Coordinate>& points, const std::vector<PointData>& data, VariogramModel model = VariogramModel::LINEAR, double a = 0, double b = 0)
		: points_(points), data_(data), model_(model), cofa_(a), cofb_(b) {}

	// ����Kriging��ֵ
	std::map<std::string, double> interpolate(const Coordinate& target);

private:
	std::vector<std::pair<double, double>> experimentalVariogram(const std::vector<Coordinate>& points, const std::vector<double>& values);
	std::pair<double, double> fitVariogramModel(const std::vector<std::pair<double, double>>& variogram);
	double variogramFunction(double h)const;
	double interpolateValue(const std::vector<double>& values, const Coordinate& target);  // �����һ��

	std::vector<Coordinate> points_;
	std::vector<PointData> data_;
	VariogramModel model_;
	double cofa_;
	double cofb_;
};
