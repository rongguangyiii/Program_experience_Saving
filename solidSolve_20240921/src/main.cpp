#include <iostream>
#include <fstream>
#include <memory>
#include "grid/include/Basegrid.h"
#include "grid/include/gridSin.h"
#include "grid/include/gridDisc.h"
#include "grid/include/gridPerturb.h"
#include "grid/include/gridUniform.h"
#include "grid/include/globalData.h"
#include "solve/include/SolidSolver.h"
#include "solve/include/SolidSolverLSG.h"
#include "solve/include/SolidSolverChain.h"
#include "solve/include/SolidSolverLSGNoDEER.h"

//outGrid(cur, name[gridIndex]);
static void outGrid(const std::shared_ptr<GridBase>& it, const std::string& name) {
	std::string Tecplotname = "../../Grid/" + name + "_Tecplot.dat";
	std::ofstream out(Tecplotname);
	it->outputToTecplot(out, Tecplotname);
	out.close();
}
//outGrid(gridVec);
static void outGrid(const std::vector<std::shared_ptr<GridBase>>& gridVec) {
	GlobalData::createFolder("../../Grid");
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
		grid = std::make_shared<GridUniform>(0.0, 1.0, 0.0, 1.0, 21, 21);//����
	}
	else if (i == 1) {
		grid = std::make_shared<GridBoundS>(0.0, 1.0, 0.0, 1.0, 21, 21, 2, 2);//����
	}
	else if (i == 2) {
		grid = std::make_shared<DiscGrid>(0.0, 1.0, 0.0, 1.0, 51, 51);//Բ��
	}
	else if (i == 3) {
		grid = std::make_shared<GridPerturb>(0.0, 1.0, 0.0, 1.0, 21, 21);//�Ŷ�
	}
	else {
		std::cout << "ERROR!\n";
	}
	return grid;
}

int main() {

	const std::shared_ptr<GridBase>& cur = genGrid(3);//uni��0, sin��1, dic��2, perturb��3
	//outGrid(cur, "perturb");
	size_t nx = cur->getXnum();
	size_t ny = cur->getYnum();
	//SolidSolver ctrlSolve(cur, nx, ny);
	//SolidSolverChain ctrlSolve(cur, nx, ny);
	SolidSolverLSG ctrlSolve(cur, nx, ny);
	//SolidSolverLSGNoDEER ctrlSolve(cur, nx, ny);
	ctrlSolve.solve();
	return 0;
}
