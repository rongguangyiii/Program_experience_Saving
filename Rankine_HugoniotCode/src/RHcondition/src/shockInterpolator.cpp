#include "RHcondition/include/shockInterpolator.h"
#include <iostream>
#include <map>
#include <string>

// 计算IDW插值-----------------------------------------------------------------------------------------
std::map<std::string, double> IDWInterpolator::interpolate(const Coordinate& targetP) const
{
	std::map<std::string, double> result;
	double numeratorTemp = 0.0, numeratorDensity = 0.0, numeratorPressure = 0.0;
	double denominator = 0.0;

	for (size_t i = 0; i < coordinates_.size(); ++i) {
		double dist = distance(targetP, coordinates_[i]);
		if (std::fabs(dist - 0) < 1E-8) {
			result["temperature"] = pointData_[i].temperature_;
			result["density"] = pointData_[i].density_;
			result["pressure"] = pointData_[i].pressure_;
			return result;  // 如果查询点恰好是已知点，直接返回已知点的属性
		}
		double weight = 1.0 / pow(dist, power_);
		numeratorTemp += weight * pointData_[i].temperature_;
		numeratorDensity += weight * pointData_[i].density_;
		numeratorPressure += weight * pointData_[i].pressure_;
		denominator += weight;
	}

	result["temperature"] = numeratorTemp / denominator;
	result["density"] = numeratorDensity / denominator;
	result["pressure"] = numeratorPressure / denominator;

	return result;
}

// 计算Kriging插值-----------------------------------------------------------------------------------------
std::map<std::string, double> KrigingInterpolator::interpolate(const Coordinate& target) {
	std::map<std::string, double> result;

	// 插值温度
	std::vector<double> temperatures;
	for (const auto& d : data_) {
		temperatures.push_back(d.temperature_);
	}
	double temp = interpolateValue(temperatures, target);
	result["temperature"] = temp;

	// 插值密度
	std::vector<double> densities;
	for (const auto& d : data_) {
		densities.push_back(d.density_);
	}
	double density = interpolateValue(densities, target);
	result["density"] = density;

	// 插值压力
	std::vector<double> pressures;
	for (const auto& d : data_) {
		pressures.push_back(d.pressure_);
	}
	double pressure = interpolateValue(pressures, target);
	result["pressure"] = pressure;

	return result;
}

double KrigingInterpolator::interpolateValue(const std::vector<double>& values, const Coordinate& target) {
	// 计算实验半变异函数
	std::vector<std::pair<double, double>> variogram = experimentalVariogram(points_, values);

	// 拟合变异函数模型参数
	std::pair<double, double> params = fitVariogramModel(variogram);
	cofa_ = params.first;
	cofb_ = params.second;

	int n = points_.size();

	// 距离矩阵
	Eigen::MatrixXd distances(n, n);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			distances(i, j) = distance(points_[i], points_[j]);
		}
	}

	// 半变异矩阵
	Eigen::MatrixXd gamma(n + 1, n + 1);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			gamma(i, j) = variogramFunction(distances(i, j));
		}
	}
	gamma.block(n, 0, 1, n).setOnes();
	gamma.block(0, n, n, 1).setOnes();
	gamma(n, n) = 0;

	// 目标点与已知点的距离
	Eigen::VectorXd gamma_vector(n + 1);
	for (int i = 0; i < n; ++i) {
		gamma_vector(i) = variogramFunction(distance(points_[i], target));
	}
	gamma_vector(n) = 1;

	// 求解权重
	Eigen::VectorXd weights = gamma.fullPivLu().solve(gamma_vector);

	// 计算插值结果
	double value = 0;
	for (int i = 0; i < n; ++i) {
		value += weights(i) * values[i];
	}

	return value;
}

// 计算实验半变异函数
std::vector<std::pair<double, double>> KrigingInterpolator::experimentalVariogram(const std::vector<Coordinate>& points, const std::vector<double>& values) {
	int n = points.size();
	std::vector<std::pair<double, double>> variogram;

	for (int i = 0; i < n; ++i) {
		for (int j = i + 1; j < n; ++j) {
			double h = distance(points[i], points[j]);
			double gamma_h = 0.5 * pow(values[i] - values[j], 2);
			variogram.push_back({ h, gamma_h });
		}
	}

	return variogram;
}

// 拟合变异函数模型参数
std::pair<double, double> KrigingInterpolator::fitVariogramModel(const std::vector<std::pair<double, double>>& variogram) {
	int n = variogram.size();
	Eigen::MatrixXd X(n, 2);
	Eigen::VectorXd y(n);

	for (int i = 0; i < n; ++i) {
		X(i, 0) = 1;
		X(i, 1) = variogram[i].first;
		y(i) = variogram[i].second;
	}

	Eigen::VectorXd beta = (X.transpose() * X).ldlt().solve(X.transpose() * y);
	double a = beta(0);
	double b = beta(1);

	return { a, b };
}

// 计算变异函数值
double KrigingInterpolator::variogramFunction(double h)const {
	switch (model_) {
	case VariogramModel::LINEAR:
		return cofa_ + cofb_ * h;
	case VariogramModel::SPHERICAL:
		if (h <= cofb_) {
			return cofa_ * (1.5 * (h / cofb_) - 0.5 * pow(h / cofb_, 3));
		}
		else {
			return cofa_;
		}
	case VariogramModel::EXPONENTIAL:
		return cofa_ * (1 - exp(-h / cofb_));
	case VariogramModel::GAUSSIAN:
		return cofa_ * (1 - exp(-pow(h / cofb_, 2)));
	default:
		return cofa_ + cofb_ * h;
	}
}

