#pragma once
#include "grid/include/baseGrid.h"

namespace Tools 
{
	void createFolder(const std::string& basePath);
	void outGrid(const std::shared_ptr<GridBase>& it, const std::string& name);
	void outGrid(const std::vector<std::shared_ptr<GridBase>>& gridVec);
	std::shared_ptr<GridBase> genGrid(size_t i);

}