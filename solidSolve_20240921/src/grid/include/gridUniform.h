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
	//���ز�����=
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
 USAGE: �ṩ���η�Χ����������������
 -------------------------------------------------------------------------------
 int main() {
	// ��һ����������
	UniformGrid grid1(0.0, 4.0, 0.0, 4.0, 21, 21);
	// �ڶ�����������
	UniformGrid grid2(6.0, 10.0, 0.0, 4.0, 21, 21);
	// ������ļ�
	std::ofstream out("uniform_zones_with_elements.dat");
	// �����һ������
	grid1.outputToTecplot("Zone 1", out);
	// ����ڶ�������
	grid2.outputToTecplot("Zone 2", out);
	out.close();
	std::cout << "Points and elements output to uniform_points_zones_with_elements.dat\n";
	return 0;
 }
 -------------------------------------------------------------------------------
 ʾ��2��

int main() {
	// ����һ�� 3x3 ������
	UniformGrid grid(0.0, 2.0, 0.0, 2.0, 3, 3);

	// ��ȡ�ھ���Ϣ
	auto neighbors = grid.getNeighbors();

	// ��ӡÿ������ھ���Ϣ
	for (const auto& [point, neighbor_points] : neighbors) {
		std::cout << "Point (" << point.x << ", " << point.y << ") has neighbors:\n";
		for (const auto& neighbor : neighbor_points) {
			std::cout << "\t(" << neighbor.x << ", " << neighbor.y <<")\n";
		}
	}

	// �������Tecplot�ļ�
	std::ofstream outfile("uniform_grid_output.dat");
	grid.outputToTecplot(outfile, "Uniform Grid");
	outfile.close();

	return 0;
}
--------------------------------------------------------------*/
