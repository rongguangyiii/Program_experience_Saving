#pragma once
#include <string>

namespace GlobalData
{
	const std::string meshType = "Uniform";//网格类型 Uniform; Semicircular;
	const std::string calmethod = "Common";//计算RHS方法 Deer(冻结系数); Common(使用计算点度量系数); Basic(无通量分裂)
	const size_t nx = 101;
	const size_t ny = 51;
	const double cfl = 0.9;      // CFL数
	const double ctrltime = 100;      // 总时间
	const size_t ctrlstep = 20000;      // 时间步数
	const double tolerance = 1e-15;  // 收敛阈值
	const size_t ctrlout = 100; //输出控制
	const size_t movestep = 1;
	const double PI = 3.141592653589793238462643383279502884197169399;//圆周率pi
	void createFolder(const std::string& basePath);
	void outGrid();
}