#include <iostream>
#include <fstream>
#include <memory>
#include "tools/include/tools.h"
#include "grid/include/Basegrid.h"
#include "grid/include/gridSin.h"
#include "grid/include/gridDisc.h"
#include "grid/include/gridPerturb.h"
#include "grid/include/gridUniform.h"
#include "grid/include/globalData.h"
#include "solve/include/flowSolveND1UP.h"
#include "solve/include/flowSolveND2UP.h"

//outGrid(cur, name[gridIndex]);
static void outGrid(const std::shared_ptr<GridBase>& it, const std::string& name) {
	std::string Tecplotname = "../../Grid/" + name + "_Tecplot.dat";
	std::ofstream out(Tecplotname);
	it->outputToTecplot(out, Tecplotname);
	out.close();
}
//outGrid(gridVec);
static void outGrid(const std::vector<std::shared_ptr<GridBase>>& gridVec) {
	Tools::createFolder("../../Grid");
	size_t num1 = 0;
	std::vector<std::string>name{ "mesh_uni", "mesh_sin" ,"mesh_dic" };
	for (const auto& it : gridVec) {
		outGrid(it, name[num1++]);
	}
}
//gen per grid
static std::shared_ptr<GridBase> genGrid(size_t i)
{
	std::shared_ptr<GridBase> grid = nullptr;
	if (i == 0) {
		grid = std::make_shared<GridUniform>(0.0, 1.0, 0.0, 1.0, 21, 21);//均匀
	}
	else if (i == 1) {
		grid = std::make_shared<GridBoundS>(0.0, 1.0, 0.0, 1.0, 21, 21, 2, 2);//正弦
		//grid = std::make_shared<GridBoundS>(0.0, 1.0, 0.0, 1.0, 41, 41, 4, 4);//正弦
		//grid = std::make_shared<GridBoundS>(0.0, 1.0, 0.0, 1.0, 81, 81, 8, 8);//正弦
		//grid = std::make_shared<GridBoundS>(0.0, 1.0, 0.0, 1.0, 161, 161, 16, 16);//正弦
	}
	else if (i == 2) {
		grid = std::make_shared<DiscGrid>(0.0, 1.0, 0.0, 1.0, 51, 51);//圆盘
	}
	else if (i == 3) {
		grid = std::make_shared<GridPerturb>(0.0, 1.0, 0.0, 1.0, 21, 21);//扰动
	}
	else {
		std::cout << "ERROR!\n";
	}
	return grid;
}

int main() {
	const std::shared_ptr<GridBase>& cur = genGrid(0);//uni：0, sin：1, dic：2, perturb：3
	//outGrid(cur, "perturb");
	size_t nx = cur->getXnum();
	size_t ny = cur->getYnum();
	//FlowSolveND1UP ctrlSolve(cur, nx, ny);
	FlowSolveND2UP ctrlSolve(cur, nx, ny);
	ctrlSolve.solve();
	return 0;
}
