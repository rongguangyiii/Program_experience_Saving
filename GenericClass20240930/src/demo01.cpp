
#include "triangle/include/triBase.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <map>

int main() {
	//TriBase test01("ring_ring");
	//test01.boostMesh();
	//--------------------------------------------------------
	TriBase test;
	test.setConstrainType("ring_ring");
	std::vector<std::vector<double>> poly_1_outter = { {1,2}, {3,2}, {7,2}, {7,4}, {8,7}, {5,8}, {2,6}, {3,4} };
	std::vector<std::vector<double>> poly_2_inner = { {4,4}, {6,3}, {6,6}, {4,6}, {5,5} };
	size_t index = 0;
	std::map<size_t, std::vector<double>>outterMap;
	std::map<size_t, std::vector<double>>innerMap;
	for (size_t i = 0; i < poly_1_outter.size(); i++){
		index++;
		outterMap[index] = poly_1_outter[i];
	}
	for (size_t i = 0; i < poly_2_inner.size(); i++){
		index++;
		innerMap[index] = poly_2_inner[i];
	}
	test.insertCtrlOutterPoly(outterMap);
	test.insertCtrlInnerPoly(innerMap);
	test.boostMesh();

	return 0;
}
