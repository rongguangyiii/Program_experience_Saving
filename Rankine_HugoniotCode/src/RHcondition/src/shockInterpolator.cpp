#include "RHcondition/include/shockInterpolator.h"
#include <iostream>
#include <map>
#include <string>

// ����IDW��ֵ-----------------------------------------------------------------------------------------
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
			return result;  // �����ѯ��ǡ������֪�㣬ֱ�ӷ�����֪�������
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

// ����Kriging��ֵ-----------------------------------------------------------------------------------------
std::map<std::string, double> KrigingInterpolator::interpolate(const Coordinate& target) {
	std::map<std::string, double> result;

	// ��ֵ�¶�
	std::vector<double> temperatures;
	for (const auto& d : data_) {
		temperatures.push_back(d.temperature_);
	}
	double temp = interpolateValue(temperatures, target);
	result["temperature"] = temp;

	// ��ֵ�ܶ�
	std::vector<double> densities;
	for (const auto& d : data_) {
		densities.push_back(d.density_);
	}
	double density = interpolateValue(densities, target);
	result["density"] = density;

	// ��ֵѹ��
	std::vector<double> pressures;
	for (const auto& d : data_) {
		pressures.push_back(d.pressure_);
	}
	double pressure = interpolateValue(pressures, target);
	result["pressure"] = pressure;

	return result;
}

double KrigingInterpolator::interpolateValue(const std::vector<double>& values, const Coordinate& target) {
	// ����ʵ�����캯��
	std::vector<std::pair<double, double>> variogram = experimentalVariogram(points_, values);

	// ��ϱ��캯��ģ�Ͳ���
	std::pair<double, double> params = fitVariogramModel(variogram);
	cofa_ = params.first;
	cofb_ = params.second;

	int n = points_.size();

	// �������
	Eigen::MatrixXd distances(n, n);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			distances(i, j) = distance(points_[i], points_[j]);
		}
	}

	// ��������
	Eigen::MatrixXd gamma(n + 1, n + 1);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			gamma(i, j) = variogramFunction(distances(i, j));
		}
	}
	gamma.block(n, 0, 1, n).setOnes();
	gamma.block(0, n, n, 1).setOnes();
	gamma(n, n) = 0;

	// Ŀ�������֪��ľ���
	Eigen::VectorXd gamma_vector(n + 1);
	for (int i = 0; i < n; ++i) {
		gamma_vector(i) = variogramFunction(distance(points_[i], target));
	}
	gamma_vector(n) = 1;

	// ���Ȩ��
	Eigen::VectorXd weights = gamma.fullPivLu().solve(gamma_vector);

	// �����ֵ���
	double value = 0;
	for (int i = 0; i < n; ++i) {
		value += weights(i) * values[i];
	}

	return value;
}

// ����ʵ�����캯��
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

// ��ϱ��캯��ģ�Ͳ���
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

// ������캯��ֵ
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

