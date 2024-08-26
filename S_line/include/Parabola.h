#include <eigen3/Eigen/Dense>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include "uniformGrid.h"

using namespace Eigen;
using namespace std;

class Parabola {
public:
	// 构造函数：左开口抛物线
	Parabola(const vector<Vector2d>& points, double step, std::string which);

	// 生成抛物线上的点
    void generatePoints();
	void outputToTecplot(ofstream& out, const string& name) const;
    vector<Coord> getPoints() const { return data_points_; }
private:
	std::string tag_;
	double a_, b_, c_, step_;
	vector<Vector2d> points_;
    vector<Coord> data_points_;


	// 计算左开口抛物线的系数
	void calLeftParabolaCoefficients(const vector<Vector2d>& points);
	// 计算右开口抛物线的系数
	void calRightParabolaCoefficients(const vector<Vector2d>& points);
};

/*-----------------------------------------------------------------------------------------------------
USGAE: 该类的用法
        给定三个点，用于生成过三个点的抛物线，目前只有向右和向左开口的抛物线
 ------------------------------------------------------------------------------------------------------
 int main() {
    // 定义抛物线的点
    vector<Vector2d> points_up = { Vector2d(0.5, 1.5), Vector2d(1, 2), Vector2d(0.5, 2.5) };
    vector<Vector2d> points_down = { Vector2d(0.5, 0.5), Vector2d(0, 1), Vector2d(0.5, 1.5) };

    // 打开一个文件用于输出到 Tecplot
    ofstream out("line_out.dat");
    out << "TITLE = \"S_LINE\"\n";
    out << "VARIABLES = \"X\", \"Y\"\n";

    // 输出开口向左的抛物线
    Parabola left_parabola(points_up, "left");
    left_parabola.outputToTecplot(out, "Left Parabola");

    // 输出开口向右的抛物线
    Parabola right_parabola(points_down, "right");
    right_parabola.outputToTecplot(out, "Right Parabola");

    out.close();
    cout << "Parabolas output to parabolas.dat\n";

    return 0;
}
-----------------------------------------------------------------------------------------------------*/