#pragma once
#include <random>
#include <vector>
#include "Coord.h"
#include "basegrid.h"
#include "coefTrans.h"


class Perturbation {
private:
    double dx_; // x方向的距离步长
    double dy_; // y方向的距离步长

    std::random_device rd_;
    std::mt19937 gen_;
    std::uniform_real_distribution<double> dist_;
public:
    // 构造函数，初始化速度a，时间步长dt，以及x和y方向的距离步长dx, dy
    Perturbation(double dx, double dy)
        : dx_(dx), dy_(dy), gen_(rd_()), dist_(0.0, dx) {}

    // 计算新坐标，使用随机速度分量ax, ay
    std::pair<double, double> calPerturbedCoord(double x0, double y0) {
        double ax, ay;
        do {
            ax = dist_(gen_);
            ay = dist_(gen_);
        } while (ax >= 0.5 * dx_ || ay >= 0.5 * dy_);

        double x1 = x0 + ax;
        double y1 = y0 + ay;
        return { x1, y1 };
    }

};

class GridPerturb : public GridBase {
public:
	GridPerturb(double x_start, double x_end, double y_start, double y_end, int x_points, int y_points)
		: x_start_(x_start), x_end_(x_end), y_start_(y_start), y_end_(y_end), level_(2),
		x_points_(x_points), y_points_(y_points), x(x_points_, std::vector<double>(y_points_, 0.0)),
		y(x_points_, std::vector<double>(y_points_, 0.0))
	{
		generatePoints();
		generateElements();
		genNeiborNode();
		coefTrans();
	}

	~GridPerturb() override = default;
	void generatePoints() override;
	void generateElements() override;
	void genNeiborNode()override;
	void coefTrans()override;
	void generatePoints_1();

	//重载操作符=
	std::shared_ptr<GridBase> operator=(const std::shared_ptr<GridBase>& other) override;
	void outputToTecplot(std::ofstream& out, const std::string& zone_title) const override;
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
