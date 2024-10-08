/*---------------------------------------------------------------------------
	ZaRan	-	A Totally Automatic CFD Software
	Copyright (C) ,Since 2020
-------------------------------------------------------------------------------
License
	This file is part of ZaRan.

!	@file	tri info
!	@brief	Generate ZaRan Grid.
!	@author	Liu Guangying.
!   @date  2024.04.20
!   @location ShenZhen
\*---------------------------------------------------------------------------*/

#include "triangle/include/cgalTypes.h"
#include "triangle/include/triEle.h"
#include "triangle/include/triBase.h"
//#include "utility/include/tools.h"
////#include "utility/include/log.h"
////#include "Basic/include/GlobalData.h"
#include <memory>
#include <vector>
#include "gridGenerate/include/coord.h"
#include "flowSolve/include/globalData.h"




/*---------------------------------------------------------------------------------------------
*  Function: size_t operator()(const TriEdge& k) const
*  Purpose: 
*      根据输入的键(TriEdge),通过FNV-1a哈希算法计算得到一个哈希值(唯一整数)
*  Description: 
*     std::unordered_map 会自动使用EdgeHash的operator()来计算edge的哈希值，并据此
*     确定在哈希表中应该存储这个键值对的位置。你不需要（也不应该）直接调用EdgeHash的operator()。
*     这是std::unordered_map的内部实现细节，它负责确保在需要时正确地使用哈希函数。
*  Usage: 
*		std::unordered_map<TriEdge, std::vector<TriEle>, EdgeHash>
*  Notice : 
*		需要重载 键(TriEdge) 对应的operator==()函数，和hash值实现函数，也就是本函数。
*  Arguments:
*    k - 输入的键，即边，由两个顶点组成
*  Returns:
*   size_t - 返回一个哈希值
---------------------------------------------------------------------------------------------*/
std::size_t EdgeHash::operator()(const std::shared_ptr<TriEdge>& k) const {
	// 使用FNV-1a哈希算法的变体，它对于整数输入通常表现良好  
	std::size_t hash = 2166136261; // FNV offset basis  
	hash ^= k->vA;
	hash *= 16777619; // FNV prime  
	hash ^= k->vB;
	hash *= 16777619; // 再次混合  
	return hash;
}

void TriBase::addEdgeTableItem(const triElePtr_& triEle_ptr)
{
	// 对每个边，将其添加到边表中
	for (const auto& edge : triEle_ptr->getTriEdgesRef()) 
	{
		if (edge->vA == edge->vB){
			//spdlog::error("当前边的TriNode出错：vA={},vB={}", edge->vA, edge->vB);
			throw std::runtime_error("Caught an exception!");
		}
		edgeTable_[edge].push_back(triEle_ptr);
	}

}

TriBase::TriBase(std::vector<CoordPtr>& outterNode, std::vector<CoordPtr>& innerNode,
	const std::string& typeName) :dividType_(typeName)
{
	if(dividType_=="ring_ring")
	{
		GlobalData::sortPointsCounterclockwise(outterNode);
		GlobalData::sortPointsCounterclockwise(innerNode);
	}
	else if (dividType_ == "ring_line")
	{
	}
	else if (dividType_ == "ring_points")
	{
		GlobalData::sortPointsCounterclockwise(innerNode);
	}

	//生成TriNodeMap，每个TriNode必须在实例化的时候就给定点类型,每个点的索引必须唯一；
	size_t index = 0;
	for (const auto& per_it : outterNode)
	{
		index++;
		const auto& currentcood = per_it;
		//const auto& index = per_it->id_;
		TriNode Point_cdt{ *currentcood, index, TriNodeType::outter };
		TriNodeMap_[index] = Point_cdt;
		cgalPointsIndexVec_.emplace_back(index);
	}
	for (const auto& per_it : innerNode)
	{
		index++;
		const auto& currentcood = per_it;
		//const auto& index = per_it->id_;
		TriNode Point_cdt{ *currentcood, index, TriNodeType::inner };
		TriNodeMap_[index] = Point_cdt;
		cgalPointsIndexVec_.emplace_back(index);
	}

}

/*----------------------------------------------------------------------------------------------
	本部分是一个简单的内凹外凹多边形；用于DEBUG
------------------------------------------------------------------------------------------------*/
TriBase::TriBase():dividType_("ring_ring")
{
	std::vector<Coord> poly_1_outter = { {1,2}, {3,2}, {7,2}, {7,4}, {8,7}, {5,8}, {2,6}, {3,4} };
	std::vector<Coord> poly_2_inner = { {4,4}, {6,3}, {6,6}, {4,6}, {5,5} };
	size_t index = 1;
	for (const auto& per_it : poly_1_outter)
	{
		TriNode Point_cdt{ per_it, index, TriNodeType::outter };
		TriNodeMap_[index] = Point_cdt;
		cgalPointsIndexVec_.emplace_back(index);
		index++;
	}
	for (const auto& per_it : poly_2_inner)
	{
		TriNode Point_cdt{ per_it, index, TriNodeType::inner };
		TriNodeMap_[index] = Point_cdt;
		cgalPointsIndexVec_.emplace_back(index);
		index++;
	}
}

/*-------------------------------------------------------------------
*  Function: updateEdgeTable()
*  Purpose: 更新边表(哈希值表)edgeTable_，第一次没有边表的话就创建哈希值表
*  Arguments:
*    null
*  Returns:
*    null
-------------------------------------------------------------------*/
void TriBase::genEdgeTable()
{
	edgeTable_.clear();
	for (auto& triEle : triEleVec_)
		addEdgeTableItem(triEle.second);
	
	return;
}

/*-------------------------------------------------------------------
*  Function: GenCompleteCell()
*  Purpose: 生成cell_,用于存储三角形和四边形的顶点索引
*  Arguments:
*    null
*  Returns:
*    null
-------------------------------------------------------------------*/
void TriBase::GenCompleteCell()
{
	for (auto& triEle : triEleVec_)
	{
		if (TriEleTag::used == triEle.second->getTriEleTag())
			continue;
		cell_.push_back(triEle.second->GetTriEleVertex());
	}

	return;
}

/*-------------------------------------------------------------------------------------------------------------
*  Function: OutputConnectGrid2Tecplot()
*  Purpose: 输出三角形和四边形网格到 Tecplot 文件中，用于查看结果。
*  Arguments:
*    pointIndexVec - 所有的点的索引。
*    cellVec - 所有的三角形和思辨，用于输出Tecplot单元元素连接表。
*    vertex2tecplotIndex - 索引映射map,用于转换节点索引从0开始。
*    title - Tecplot文件名。
*  Usage:
*		OutputConnectGrid2Tecplot(vertex2tecplotIndex_.first, cellVec, vertex2tecplotIndex_.second, "Tri_ini");
*  Returns:
*    null
--------------------------------------------------------------------------------------------------------------*/
void TriBase::OutputConnectGrid2Tecplot(const std::vector<size_t> pointIndexVec, cellVec_& cellVec,
	const std::unordered_map<size_t, size_t>& vertex2tecplotIndex, const std::string title)const
{
	std::string outfilename = "Data_" + title + "_grid" + ".dat";
	std::string basePath = "tempresult";
	GlobalData::createFolder(basePath);
	// 创建一个 Tecplot 文件

	std::ofstream tecplotFile(basePath + "/" + outfilename);
	tecplotFile << "TITLE = \"Triangle Output   LGYYYY\"\n";
	tecplotFile << "VARIABLES = \"X\", \"Y\"\n";

	// 写入顶点坐标
	tecplotFile << "ZONE T=\"Triangles\", N=" << pointIndexVec.size()
		<< ", E=" << cellVec.size() << ", F=FEPOINT, ET=QUADRILATERAL\n";
	for (const auto& vit : pointIndexVec)
	{
		const auto& p = index2Coor(vit);
		tecplotFile << p.x << " " << p.y << "\n";
		//spdlog::debug("输出节点坐标：{},{}", p.x(), p.y());
	}

	// 写入单元信息
	for (const auto& cell : cellVec)
	{
		if (cell.size() == 3) { // 三角形单元
			tecplotFile << vertex2tecplotIndex.at(cell[0]) + 1 << " " << vertex2tecplotIndex.at(cell[1]) + 1 << " "
				<< vertex2tecplotIndex.at(cell[2]) + 1 << " " << vertex2tecplotIndex.at(cell[0]) + 1 << std::endl;
		}
		else if (cell.size() == 4) { // 四边形单元
			tecplotFile << vertex2tecplotIndex.at(cell[0]) + 1 << " " << vertex2tecplotIndex.at(cell[1]) + 1 << " "
				<< vertex2tecplotIndex.at(cell[2]) + 1 << " " << vertex2tecplotIndex.at(cell[3]) + 1 << std::endl;
		}
		else {
			//spdlog::error("函数：OutputConnectGrid2Tecplot(): 单元的顶点数量不对！！！");
		}
	}

	tecplotFile.close();

	return;
}

/*-------------------------------------------------------------------
*  Function: Gen_initialCDT
*  Purpose: 生成初始三角形网格并输出到Tecplot文件中，用于查看结果。
*  Arguments:
*    null
*  Returns:
*    null
-------------------------------------------------------------------*/
void TriBase::genCDT()
{
	//step 1: 设置CDT约束边界和三角形筛选边界。
	std::vector<Point_2> cgal_outerPointsVec;
	std::vector<Point_2> cgal_innerPointsVec;
	std::vector<Coord> poly_outerCoor;
	std::vector<Coord> poly_innerCoor;
	//TriNodeMap_是没有顺序的，在TriNodeMap_中循环就会导致控制边界也没有顺序了，插入边界时候出错。
	//在cgalPointsIndexVec_中循环就保证了顺序。
	for (const auto& perit: cgalPointsIndexVec_)
	{
		const auto& cur_trinode = TriNodeMap_.at(perit);
		if (cur_trinode.gettype() == TriNodeType::outter)
		{
			cgal_outerPointsVec.push_back({ cur_trinode.getcoord().x ,cur_trinode.getcoord().y });
			poly_outerCoor.emplace_back(cur_trinode.getcoord());
		}
		else if (cur_trinode.gettype() == TriNodeType::inner)
		{
			cgal_innerPointsVec.push_back({ cur_trinode.getcoord().x ,cur_trinode.getcoord().y });
			poly_innerCoor.emplace_back(cur_trinode.getcoord());
		}
		else
		{
			//spdlog::error("current tri node type is error!!!");
			throw std::runtime_error("Caught an exception!");
		}
	}

	//step 2: 合并两个点集  cgalPointsIndexVec_的顺序决定了cgal_PointsVec的顺序！！！
	std::vector<Point_2> cgal_PointsVec(cgal_outerPointsVec.begin(), cgal_outerPointsVec.end());
	cgal_PointsVec.insert(cgal_PointsVec.end(), cgal_innerPointsVec.begin(), cgal_innerPointsVec.end());

	//step 3: 创建一个约束三角剖分
	CDT cdt;

	//step 4: 将点和索引关联起来
	std::unordered_map<size_t, size_t> vertex2tecplotIndex;//将节点索引映射到局部从0开始的索引，用于输出Tecplot文件
	std::map<Vertex_handle, size_t> vertex_to_index;//将CDT索引映射到节点索引
	for (size_t iPoint = 0; iPoint != cgalPointsIndexVec_.size(); ++iPoint)
	{
		Vertex_handle vh = cdt.insert(cgal_PointsVec.at(iPoint));
		vertex_to_index[vh] = cgalPointsIndexVec_.at(iPoint);
		vertex2tecplotIndex[cgalPointsIndexVec_.at(iPoint)] = iPoint;
	}
	vertex2tecplotIndex_.first = cgalPointsIndexVec_;//节点索引
	vertex2tecplotIndex_.second = vertex2tecplotIndex;//局部索引
	
	/*-----------------------------------------------------------------
	 * Note：CGAL库的insert_constraint函数，用于插入约束边界；
	 *	     通常需要两个类型为Point_2的参数。表示插入约束线段的两个端点，
	 *       当有多个线段插入时，选段首位相连的节点必须是有序的。
	 *-----------------------------------------------------------------*/
	//step 5: 插入约束边界,因为是两个圆环约所以束边界要分别插入！！！
	if (dividType_ == "ring_ring")
	{
		for (size_t i = 0; i < cgal_outerPointsVec.size(); ++i)
			cdt.insert_constraint(cgal_outerPointsVec[i], cgal_outerPointsVec[(i + 1) % cgal_outerPointsVec.size()]);
		for (size_t i = 0; i < cgal_innerPointsVec.size(); ++i)
			cdt.insert_constraint(cgal_innerPointsVec[i], cgal_innerPointsVec[(i + 1) % cgal_innerPointsVec.size()]);
	}
	else if (dividType_ == "ring_line")
	{
		//for (size_t i = 0; i < cgal_outerPointsVec.size(); ++i)
		//	cdt.insert_constraint(cgal_outerPointsVec[i], cgal_outerPointsVec[(i + 1) % cgal_outerPointsVec.size()]);
		for (size_t i = 0; i < cgal_innerPointsVec.size()-1; ++i)
			cdt.insert_constraint(cgal_innerPointsVec[i], cgal_innerPointsVec[i + 1]);
	}
	else if (dividType_ == "ring_points")
	{
		for (size_t i = 0; i < cgal_innerPointsVec.size(); ++i)
			cdt.insert_constraint(cgal_innerPointsVec[i], cgal_innerPointsVec[(i + 1) % cgal_innerPointsVec.size()]);
	}

	//step 6: 写入三角形索引
	size_t triIndex = 0;
	for (Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit)
	{
		size_t v0 = vertex_to_index[fit->vertex(0)];
		size_t v1 = vertex_to_index[fit->vertex(1)];
		size_t v2 = vertex_to_index[fit->vertex(2)];
		Coord p0 = index2Coor(v0);
		Coord p1 = index2Coor(v1);
		Coord p2 = index2Coor(v2);
		//step 1 : 判断三角形的三个顶点是否在一条直线上,
		if (Tools::isCollinear(p0, p1, p2))
			continue;
		Coord triCenter = (p0 += p1 += p2) /= 3.0;
		if (dividType_ == "ring_ring")
		{
			//step 3 : 判断三角形是否在<外层>边界多边形<外部>，如果<在外部>则删除
			if (!Tools::isPointInPolygon_moreOuter(triCenter, poly_outerCoor))
				continue;
			//step 2 : 判断三角形是否在<内层>边界多边形<内部>，如果<在内部>则删除
			if (Tools::isPointInPolygon_moreInner(triCenter, poly_innerCoor))
				continue;
		}
		else if (dividType_ == "ring_points") {
			//step 2 : 判断三角形是否在<内层>边界多边形<内部>，如果<在内部>则删除
			if (Tools::isPointInPolygon_moreInner(triCenter, poly_innerCoor))
				continue;
		}
		else if (dividType_ == "ring_line") {
			//step 3 : 判断三角形是否在<外层>边界多边形<外部>，如果<在外部>则删除
			if (!Tools::isPointInPolygon_moreOuter(triCenter, poly_outerCoor))
				continue;
		}
		std::vector<TriNode> triNodes;
		triNodes.emplace_back(TriNodeMap_.at(v0));
		triNodes.emplace_back(TriNodeMap_.at(v1));
		triNodes.emplace_back(TriNodeMap_.at(v2));
		TriEle triEle(triNodes, triIndex++);
		triEle.setTriCenter(triCenter);
		std::shared_ptr<TriEle> triEle_ptr = std::make_shared<TriEle>(triEle);
		addTriEle(triEle_ptr);
	}

	//step 7: 输出初始三角形网格
	std::vector< std::vector<size_t> > cellVec;
	for (const auto& perTri : triEleVec_)
		cellVec.push_back(perTri.second->GetTriEleVertex());
	OutputConnectGrid2Tecplot(vertex2tecplotIndex_.first, cellVec, vertex2tecplotIndex_.second, "Tri_ini");

	//step ... 接下来可以处理三角形信息了

	return;
}

/*-------------------------------------------------------------------
*  Function: filterTriGrid
*  Purpose: 筛选三角形网格，删除在内部多边形以内的和外部多边形以外的三角形。
*  Arguments:
*    null
*  Returns:
*    null
-------------------------------------------------------------------*/
void TriBase::filterTri()
{
	const auto& triEleVEcPtr = GetTriEleVec();
	for (const auto& perit : triEleVEcPtr)
	{
		auto& currentTri = perit.second;
		//3个顶点都在外部边界上或内部边界上
		if (currentTri->getTriEleType() == TriEleType::Inner_only ||
			currentTri->getTriEleType() == TriEleType::Outer_only) {
			addNotIdealTri(*currentTri);
		}
	}
	return;
}

/*-------------------------------------------------------------------
*  Function: reconstructNotIdealTriGrid
*  Purpose: 对不符合要求的，但在内环和外环之间的三角形进行处理，
*			重新连接相邻三角形网格。
*  Arguments:
*    null
*  Returns:
*    null
-------------------------------------------------------------------*/
void TriBase::reconstructNotIdealTriGrid()
{
	//step 0: 筛选需要重构的三角形网格
	filterTri();
	std::vector<TriEle>& notIdealTriVec = getNotIdealTriVec();
	//step 0:如果没有不理想的三角形，直接返回
	if (0 == notIdealTriVec.size())
		return;

	//step 1: 生成边表
	genEdgeTable(); 
	//step 2: 重构非理想三角形
	std::vector<size_t> toDeleteTriIndex;
	for (size_t inot = 0; inot < notIdealTriVec.size(); inot++)
	{
		auto& currentTri = notIdealTriVec.at(inot);
		for (const auto& edge : currentTri.getTriEdgesRef())
		{
			auto it = getEdgeTable().find(edge);
			if (it->second.size() == 1)
				continue;

			const std::vector<std::shared_ptr<TriEle>>& triangles = it->second;
			//为防止误删三角形，先将要删除的三角形索引存储起来
			for(auto& var : triangles)
				toDeleteTriIndex.push_back(var->GetTriEleIndex());
			reconstructNotIdealTriGrid(*edge, triangles);
		}
	}
	//step 3: 删除部分三角形并同时更新边表
	for (size_t cc : toDeleteTriIndex)
		removeTriEle(cc);

	//step 4: 输出网格
	std::vector< std::vector<size_t> > cellVec;
	for (const auto& perTri : GetTriEleVec())
		cellVec.push_back(perTri.second->GetTriEleVertex());
	OutputConnectGrid2Tecplot(vertex2tecplotIndex_.first, cellVec, vertex2tecplotIndex_.second, "Tri_reconstruct");
	return;
}

/*-------------------------------------------------------------------
*  Function: reconstructNotIdealTriGrid(edge,triangles)
*  Purpose: 对当前共边的两个三角形重新连接。
*  Arguments:
*    edge - 共边
*   triangles - 共边的两个三角形
*  Returns:
*    null
-------------------------------------------------------------------*/
void TriBase::reconstructNotIdealTriGrid(const TriEdge& edge, const std::vector<triElePtr_>& tri_Ptr)
{
	std::vector<size_t> tri_1_vertex = tri_Ptr.at(0)->GetTriEleVertex();
	std::vector<size_t> tri_2_vertex = tri_Ptr.at(1)->GetTriEleVertex();
	size_t coindex_a = edge.vA;
	size_t coindex_b = edge.vB;
	size_t tri_1_diff = 0;
	size_t tri_2_diff = 0;
	//找到第三个和coindex_a，coindex_b不同的顶点
	for (auto it : tri_1_vertex)
	{
		if (it != coindex_a && it != coindex_b)
			tri_1_diff = it;
	}
	for (auto it : tri_2_vertex)
	{
		if (it != coindex_a && it != coindex_b)
			tri_2_diff = it;
	}
	std::vector<size_t> newTri_1{ tri_1_diff, tri_2_diff, coindex_a };
	std::vector<size_t> newTri_2{ tri_1_diff, tri_2_diff, coindex_b };
	std::vector<std::pair<Coord, size_t >> newTri_1_pair;
	std::vector<std::pair<Coord, size_t >> newTri_2_pair;
	Coord average_1;
	Coord average_2;
	for (auto it : newTri_1)
	{
		Coord temcoor = index2Coor(it);
		newTri_1_pair.push_back(std::make_pair(temcoor, it));
		average_1 += (temcoor /= newTri_1.size());
	}
	for (auto it : newTri_2)
	{
		Coord temcoor = index2Coor(it);
		newTri_2_pair.push_back(std::make_pair(temcoor, it));
		average_2 += (temcoor /= newTri_2.size());
	}
	Tools::reorderPointsIndex(newTri_1_pair, average_1);//按逆时针重新排序
	Tools::reorderPointsIndex(newTri_2_pair, average_2);
	newTri_1.clear();
	newTri_2.clear();
	for (auto& it : newTri_1_pair)
		newTri_1.push_back(it.second);
	for (auto& it : newTri_2_pair)
		newTri_2.push_back(it.second);

	std::vector<TriNode> triNodes_1;
	triNodes_1.emplace_back(TriNodeMap_.at(newTri_1.at(0)));
	triNodes_1.emplace_back(TriNodeMap_.at(newTri_1.at(1)));
	triNodes_1.emplace_back(TriNodeMap_.at(newTri_1.at(2)));
	TriEle triEle_1(triNodes_1, GetTriEleVec().size());
	triEle_1.setTriCenter(average_1);
	const auto& current_triptr = std::make_shared<TriEle>(triEle_1);
	addTriEle(current_triptr);
	//AddTriEle(triEle_1);
	addEdgeTableItem(current_triptr);

	std::vector<TriNode> triNodes_2;
	triNodes_2.emplace_back(TriNodeMap_.at(newTri_2.at(0)));
	triNodes_2.emplace_back(TriNodeMap_.at(newTri_2.at(1)));
	triNodes_2.emplace_back(TriNodeMap_.at(newTri_2.at(2)));
	TriEle triEle_2(triNodes_2, GetTriEleVec().size());
	triEle_2.setTriCenter(average_2);
	const auto& another_ptr = std::make_shared<TriEle>(triEle_2);
	addTriEle(another_ptr);
	//AddTriEle(triEle_2);
	addEdgeTableItem(another_ptr);

	return;
}

/*-------------------------------------------------------------------
*  Function: GenUnstructConnectElement_Tri2quad
*  Purpose: 生成背景网格和非结构网格的连接部分——用部分CDT网格和四边形网格填充
*  Arguments:
*    null
*  Returns:
*    null
-------------------------------------------------------------------*/
void TriBase::genTriAndQuad()
{
	//step 0 :判断是否要生成边表
	if (0 == edgeTable_.size())
		genEdgeTable();

	//step 1 :三角形转为四边形
	for (const auto& perit : edgeTable_)
	{
		auto& cur_edge = perit.first;
		auto& cur_triVec = perit.second;
		if (cur_triVec.size() == 0 || cur_triVec.size() >= 3) {
			//spdlog::error("边表出错！！！函数GenTriAndQuad() 835");
			throw std::runtime_error("Caught an exception!");
		}
		//当前边是边界边。
		if (cur_triVec.size() == 1)
			continue;
		//当前边对应的两个三角形至少有一个已经用于连接过。
		if (cur_triVec[0]->getTriEleTag() == TriEleTag::used ||
			cur_triVec[1]->getTriEleTag() == TriEleTag::used)
			continue;
		//当前边对应的两个三角形类型都是一样的。
		if (cur_triVec[0]->getTriEleType() == cur_triVec[1]->getTriEleType())
			continue;

		//如果当前三角形有符合要求的相邻三角形，就进行重构为四边形，并且同时标记相邻三角形为TriEleTag::used。
		genTriAndQuad(*cur_edge, cur_triVec);

	}//end for edge

	//step 2 :所有三角形循环结束后，再处理仍然未被使用的三角形
	GenCompleteCell();
	//step 3 :输出网格
	OutputConnectGrid2Tecplot(vertex2tecplotIndex_.first, getCell(), vertex2tecplotIndex_.second, "Tri_quad_final");

	return;
}

/*-------------------------------------------------------------------
*  Function: GenUnstructConnectElement_Tri2quad
*  Purpose: 生成背景网格和非结构网格的连接部分——用部分CDT网格和四边形网格填充
*  Arguments:
*    edge - 共边
*   triangles - 共边的两个三角形
*  Returns:
*    null
-------------------------------------------------------------------*/
void TriBase::genTriAndQuad(const TriEdge& edge, const std::vector<triElePtr_>& triangles)
{
	std::vector<size_t> tri_1_vertex = triangles.at(0)->GetTriEleVertex();
	std::vector<size_t> tri_2_vertex = triangles.at(1)->GetTriEleVertex();
	size_t coindex_a = edge.vA;
	size_t coindex_b = edge.vB;
	size_t tri_1_diff = 0;
	size_t tri_2_diff = 0;
	//找到第三个和coindex_a，coindex_b不同的顶点
	for (auto it : tri_1_vertex)
	{
		if (it != coindex_a && it != coindex_b)
			tri_1_diff = it;
	}
	for (auto it : tri_2_vertex)
	{
		if (it != coindex_a && it != coindex_b)
			tri_2_diff = it;
	}

	std::vector<size_t> quad{ tri_1_diff, tri_2_diff, coindex_a, coindex_b };
	std::vector<std::pair<Coord, size_t >> quad_pair;
	Coord average;
	for (auto it : quad)
	{
		Coord temcoor = index2Coor(it);
		quad_pair.push_back(std::make_pair(temcoor, it));
		average += (temcoor /= quad.size());
	}

	Tools::reorderPointsIndex(quad_pair, average);//按逆时针重新排序
	quad.clear();
	for (auto& it : quad_pair)
		quad.push_back(it.second);

	getCell().push_back(quad);//第一次getCell是cell_是空的，在此添加cell_

	//标记相邻三角形为TriEleTag::used
	for (auto& var : triangles)
	{
		if (var->getTriEleTag() == TriEleTag::used){
			//spdlog::error("索引为{}的三角形转换四边形出错!", var->GetTriEleIndex());
			throw std::runtime_error("Caught an exception!");
		}
		var->setTriEleTag(TriEleTag::used);
	}

	return;
}

void TriBase::gencell(const triEleMap_& content)
{
	//cell_.clear();
	for (auto& triEle : content)
		cell_.push_back(triEle.second->GetTriEleVertex());
	return;
}

void TriBase::boostMesh()
{
	genCDT();
	//reconstructNotIdealTriGrid();//:Q重构非理想三角形是必须的吗？A感觉不是待后续验证...
	//std::string method = GlobalData::GetString("ConnectMethod");
	std::string method = "HYB";//用于DEBUG
	if (method == "CDT")
		gencell(GetTriEleVec());
	else if (method == "HYB")
		genTriAndQuad();
	else
		//spdlog::error("gencell()出错：key word error！");
	return;
}

void TriBase::addTriEle(const triElePtr_& triEle)
{ 
	size_t index = triEle->GetTriEleIndex();
	triEleVec_.insert_or_assign(index, triEle); // 无论元素是否已经存在，都可以正确插入或更新。
	//triEleVec_.insert_or_assign(index, std::move(triEle)); // 无论元素是否已经存在，都可以正确插入或更新。
};

void TriBase::removeTriFromEdgeTable(const triElePtr_& triEle_ptr) {
	// 对每个边，从边表中移除triEle_ptr
	for (const auto& edge : triEle_ptr->getTriEdgesRef()) {
		auto it = edgeTable_.find(edge);
		if (it != edgeTable_.end()) {
			auto& triEleList = it->second;
			triEleList.erase(std::remove(triEleList.begin(), triEleList.end(), triEle_ptr), triEleList.end());
			if (triEleList.empty()) {
				edgeTable_.erase(it); // 如果列表为空，移除该边条目
			}
		}
	}
}

void TriBase::removeTriEle(size_t triEleKey) {
	auto it = triEleVec_.find(triEleKey);
	if (it != triEleVec_.end()) {
		auto& triEle_ptr = it->second; // 使用现有的共享指针
		removeTriFromEdgeTable(triEle_ptr); // 首先从边表中移除
		triEleVec_.erase(it); // 然后从 triEleVec_ 中移除
	}
	else {
		throw std::runtime_error("TriEle not found in triEleVec_");
	}
}