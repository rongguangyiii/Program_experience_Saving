#pragma once
#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <unordered_map>
#include <set>
#include "Coord.h"
#include "basegrid.h"
#include "coefTrans.h"

class GridUniform : public GridBase {
public:
	GridUniform(double x_start, double x_end, double y_start, double y_end, int x_points, int y_points)
		: x_start_(x_start), x_end_(x_end), y_start_(y_start), y_end_(y_end), level_(0),
		x_points_(x_points), y_points_(y_points), x(x_points_, std::vector<double>(y_points_, 0.0)),
		y(x_points_, std::vector<double>(y_points_, 0.0))
	{
		generatePoints();
		generateElements();
		genNeiborNode();
		coefTrans();
	}

	~GridUniform() override = default;
	void generatePoints() override;
	void generateElements() override;
	void genNeiborNode()override;
	void coefTrans()override;
	//重载操作符=
	std::shared_ptr<GridBase> operator=(const std::shared_ptr<GridBase>& other) override;
	void outputToTecplot(std::ofstream& out, const std::string& zone_title) const override;
	void outputToTecplot_VCFEM(std::ofstream& out, const std::string& zone_title) const override;
	size_t getXnum() const { return x_points_; }
	size_t getYnum() const { return y_points_; }
	double getXstart() const { return x_start_; }
	double getXend() const { return x_end_; }
	double getYstart() const { return y_start_; }
	double getYend() const { return y_end_; }
	std::vector<std::vector<double>> getXCoord() const { return x; }
	std::vector<std::vector<double>> getYCoord() const { return y; }
	const std::vector<CoefTrans>& getCoefTransVec() const { return coeftransVec_; }
private:
	double x_start_, x_end_, y_start_, y_end_;
	size_t x_points_, y_points_, level_;
	std::vector<std::vector<double>> x;
	std::vector<std::vector<double>> y;
	std::vector<CoefTrans> coeftransVec_;
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
	// 输出第一个区域
	grid1.outputToTecplot("Zone 1", out);
	// 输出第二个区域
	grid2.outputToTecplot("Zone 2", out);
	out.close();
	std::cout << "Points and elements output to uniform_points_zones_with_elements.dat\n";
	return 0;
 }
 -------------------------------------------------------------------------------
 示例2：

int main() {
	// 创建一个 3x3 的网格
	UniformGrid grid(0.0, 2.0, 0.0, 2.0, 3, 3);

	// 获取邻居信息
	auto neighbors = grid.getNeighbors();

	// 打印每个点的邻居信息
	for (const auto& [point, neighbor_points] : neighbors) {
		std::cout << "Point (" << point.x << ", " << point.y << ") has neighbors:\n";
		for (const auto& neighbor : neighbor_points) {
			std::cout << "\t(" << neighbor.x << ", " << neighbor.y <<")\n";
		}
	}

	// 输出网格到Tecplot文件
	std::ofstream outfile("uniform_grid_output.dat");
	grid.outputToTecplot(outfile, "Uniform Grid");
	outfile.close();

	return 0;
}
--------------------------------------------------------------*/
