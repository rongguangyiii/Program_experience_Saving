#include "tools/include/tools.h"
#include "grid/include/gridSin.h"
#include "grid/include/gridDisc.h"
#include "grid/include/gridPerturb.h"
#include "grid/include/gridUniform.h"
#include <fstream>
#include <filesystem> 

void Tools::createFolder(const std::string& basePath)
{
	std::filesystem::path folderPath = basePath;
	if (!std::filesystem::exists(folderPath))
		std::filesystem::create_directories(folderPath);
}
//outGrid(gridVec);
void Tools::outGrid(const std::vector<std::shared_ptr<GridBase>>& gridVec) {
	createFolder("../../Grid");
	size_t num1 = 0;
	std::vector<std::string>name{ "mesh_uni", "mesh_sin" ,"mesh_dic" };
	for (const auto& it : gridVec) {
		outGrid(it, name[num1++]);
	}
}
//outGrid(cur, name[gridIndex]);
void Tools::outGrid(const std::shared_ptr<GridBase>& it, const std::string& name) {
	createFolder("../../Grid");
	std::string Tecplotname = "../../Grid/" + name + "_Tecplot.dat";
	std::ofstream out(Tecplotname);
	it->outputToTecplot(out, Tecplotname);
	out.close();
}
//gen per grid
std::shared_ptr<GridBase> Tools::genGrid(size_t i)
{
	std::shared_ptr<GridBase> grid = nullptr;
	if (i == 0) {
		grid = std::make_shared<GridUniform>(0.0, 1.0, 0.0, 1.0, 21, 21);//æ˘‘»
	}
	else if (i == 1) {
		grid = std::make_shared<GridBoundS>(0.0, 1.0, 0.0, 1.0, 21, 21, 2, 2);//’˝œ“
	}
	else if (i == 2) {
		grid = std::make_shared<DiscGrid>(0.0, 1.0, 0.0, 1.0, 51, 51);//‘≤≈Ã
	}
	else if (i == 3) {
		grid = std::make_shared<GridPerturb>(0.0, 1.0, 0.0, 1.0, 21, 21);//»≈∂Ø
	}
	else {
		std::cout << "ERROR!\n";
	}
	return grid;
}