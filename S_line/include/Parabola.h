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
	// ���캯�����󿪿�������
	Parabola(const vector<Vector2d>& points, double step, std::string which);

	// �����������ϵĵ�
    void generatePoints();
	void outputToTecplot(ofstream& out, const string& name) const;
    vector<Coord> getPoints() const { return data_points_; }
private:
	std::string tag_;
	double a_, b_, c_, step_;
	vector<Vector2d> points_;
    vector<Coord> data_points_;


	// �����󿪿������ߵ�ϵ��
	void calLeftParabolaCoefficients(const vector<Vector2d>& points);
	// �����ҿ��������ߵ�ϵ��
	void calRightParabolaCoefficients(const vector<Vector2d>& points);
};

/*-----------------------------------------------------------------------------------------------------
USGAE: ������÷�
        ���������㣬�������ɹ�������������ߣ�Ŀǰֻ�����Һ����󿪿ڵ�������
 ------------------------------------------------------------------------------------------------------
 int main() {
    // ���������ߵĵ�
    vector<Vector2d> points_up = { Vector2d(0.5, 1.5), Vector2d(1, 2), Vector2d(0.5, 2.5) };
    vector<Vector2d> points_down = { Vector2d(0.5, 0.5), Vector2d(0, 1), Vector2d(0.5, 1.5) };

    // ��һ���ļ���������� Tecplot
    ofstream out("line_out.dat");
    out << "TITLE = \"S_LINE\"\n";
    out << "VARIABLES = \"X\", \"Y\"\n";

    // ������������������
    Parabola left_parabola(points_up, "left");
    left_parabola.outputToTecplot(out, "Left Parabola");

    // ����������ҵ�������
    Parabola right_parabola(points_down, "right");
    right_parabola.outputToTecplot(out, "Right Parabola");

    out.close();
    cout << "Parabolas output to parabolas.dat\n";

    return 0;
}
-----------------------------------------------------------------------------------------------------*/