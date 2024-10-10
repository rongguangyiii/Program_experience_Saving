#include <iostream>
#include <fstream>
#include <memory>
#include "grid/include/Basegrid.h"
#include "grid/include/gridSin.h"
#include "grid/include/gridDisc.h"
#include "grid/include/gridPerturb.h"
#include "grid/include/gridUniform.h"
#include "grid/include/globalData.h"
#include "Tools/include/tools.h"

int main() {

	const std::shared_ptr<GridBase>& cur = Tools::genGrid(0);//uni£º0, sin£º1, dic£º2, perturb£º3
	Tools::outGrid(cur, "uni");
	return 0;
}
