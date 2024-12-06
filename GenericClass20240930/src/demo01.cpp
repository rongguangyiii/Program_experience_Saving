
#include "triangle/include/triBase.h"
#include <vector>

int main() {
	//--------------------------------------------------------
	//为了构建正确的CDT，这里提供三种类型，outter-inner:
	//ring_ring;  point_point; point_line; point_ring;
	//完善CDT,主要增加了拓扑优化，反馈更新背景网格点的功能
	//--------------------------------------------------------
	TriBase test;
	test.setConstrainType(DataType::ring_ring);
	std::vector<TriNode>poly_1_outter = { {1,2,0,1}, {3,2,0,2}, {7,2,0,3}, {7,4,0,4}, {8,7,0,5}, {5,8,0,6}, {2,6,0,7}, {3,4,0,8} };
	std::vector<TriNode>poly_2_inner = { {4,4,0,9}, {6,3,0,10}, {6,6,0,11}, {4,6,0,12}, {5,5,0,13} };
	//test.insertPoint(poly_1_outter, TriNodeType::outter);
	//test.insertPoint(poly_2_inner, TriNodeType::inner);
	test.insertCtrlPoly(poly_1_outter,TriNodeType::outter);
	test.insertCtrlPoly(poly_2_inner, TriNodeType::inner);
	test.setfileOnOff(true);
	test.boostMesh();

	return 0;

}
