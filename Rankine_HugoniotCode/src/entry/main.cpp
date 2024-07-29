
#include "RHcondition/include/RHcondition.h"
#include "RHcondition/include/shockInterpolator.h"
#include "Coordinate/include/Coordinate.h"
#include <iostream>

//int main() {
//    double rho1 = 1.225; // 激波前密度，单位 kg/m^3
//    double u1 = 340; // 激波前速度 u 分量，单位 m/s
//    double v1 = 0; // 激波前速度 v 分量，单位 m/s
//    double p1 = 101325; // 激波前压力，单位 Pa
//    double us = 500; // 激波速度，单位 m/s
//
//    RHcondition rh(rho1, u1, v1, p1, us);
//    rh.calculatePostShockState();
//    rh.printPostShockState();
//    return 0;
//}


int main() {
	// 已知点的坐标（x, y, z）和对应的值
	std::vector<Coordinate> coords = {
		{2.0, 1.0, 0.0},
		{3.0, 2.0, 0.0},
		{4.0, 2.0, 0.0},
		{2.5, 2.5, 0.0},
		{3.5, 2.8, 0.0},
		{3.0, 3.0, 0.0},
		{4.0, 3.0, 0.0},
		{4.0, 4.0, 0.0}
	};

	std::vector<PointData> data = {
		{26.0, 3.0, 150.0},
		{18.0, 2.8, 87.0},
		{12.7, 5.1, 64.0},
		{28.0, 4.5, 120.0},
		{14.6, 7.9, 131.0},
		{22.3, 6.3, 116.0},
		{13.1, 3.2, 155.0},
		{16.0, 1.4, 74.0}
	};

	//IDW--------------------------------------------------------------------
	// 创建IDWInterpolator对象
	//IDWInterpolator interpolator(coords, data, 2.0);
	//Krig--------------------------------------------------------------------
	//创建克里金插值对象
	KrigingInterpolator interpolator(coords, data, VariogramModel::SPHERICAL);

	// 9号点的坐标
	Coordinate target(5.0, 2.6, 0.0);

	// 计算9号点的温度
	std::map<std::string, double> interpolated_data = interpolator.interpolate(target);

	std::cout << "The interpolated data at point (5.0, 2.6，0.0) is:" << std::endl;
	std::cout << "Temperature: " << interpolated_data["temperature"] << "°C" << std::endl;
	std::cout << "Density: " << interpolated_data["density"] << std::endl;
	std::cout << "Pressure: " << interpolated_data["pressure"] << "Pa" << std::endl;

	return 0;
}
