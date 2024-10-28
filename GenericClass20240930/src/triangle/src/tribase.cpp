/*---------------------------------------------------------------------------
	"A generic class that invokes the CGAL library to generate triangle connectivity information."
	Copyright (C) ,Since 2024
-------------------------------------------------------------------------------
License

!	@file	geometry_types.h
!	@brief	CGAL link typedefs.
!	@author	Liu Guangying.
!   @date  2024.09.30
!   @location DaLian
\*---------------------------------------------------------------------------*/
#include "triangle/include/cgalTypes.h"
#include "triangle/include/triEle.h"
#include "triangle/include/triBase.h"
#include <memory>
#include <vector>

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
		if (edge->vA == edge->vB) {
			//spdlog::error("当前边的TriNode出错：vA={},vB={}", edge->vA, edge->vB);
			throw std::runtime_error("Caught an exception!");
		}
		edgeTable_[edge].push_back(triEle_ptr);
	}

}

TriBase::TriBase(std::vector<TriNode>& outterNode, std::vector<TriNode>& innerNode,
	const std::string& typeName) :dividType_(typeName)
{
	if (dividType_ == "ring_ring")
	{
		reorderNodes(outterNode);
		reorderNodes(innerNode);
	}

	//生成TriNodeMap，每个TriNode必须在实例化的时候就给定点类型,每个点的索引必须唯一；
	for (const auto& per_it : outterNode)
	{
		const size_t& index = per_it.getindex();
		TriNodeMap_[index] = per_it;
		cgalPointsIndexVec_.emplace_back(index);
	}
	for (const auto& per_it : innerNode)
	{
		const size_t& index = per_it.getindex();
		TriNodeMap_[index] = per_it;
		cgalPointsIndexVec_.emplace_back(index);
	}

}

void TriBase::insertCtrlOutterPoly(std::map<size_t, std::vector<double>>CoorMap)
{
	std::vector<TriNode> nodeVec;
	for (const auto& it : CoorMap)
	{
		size_t index = it.first;
		std::vector<double> coor = it.second;
		double x = 0, y = 0, z = 0;
		if (coor.size() == 2) {
			x = coor[0], y = coor[1];
		}
		else if (coor.size() == 3) {
			x = coor[0], y = coor[1], z = coor[2];
		}
		nodeVec.emplace_back(TriNode(x, y, z, index, TriNodeType::outter));
	}
	reorderNodes(nodeVec);
	for (size_t i = 0; i < nodeVec.size(); i++)
	{
		size_t index = nodeVec[i].getindex();
		TriNodeMap_[index] = nodeVec[i];
		cgalPointsIndexVec_.emplace_back(index);
	}

}
void TriBase::insertCtrlInnerPoly(std::map<size_t, std::vector<double>>CoorMap)
{
	std::vector<TriNode> nodeVec;
	for (const auto& it : CoorMap)
	{
		size_t index = it.first;
		std::vector<double> coor = it.second;
		double x = 0, y = 0, z = 0;
		if (coor.size() == 2) {
			x = coor[0], y = coor[1];
		}
		else if (coor.size() == 3) {
			x = coor[0], y = coor[1], z = coor[2];
		}
		nodeVec.emplace_back(TriNode(x, y, z, index, TriNodeType::inner));
	}
	reorderNodes(nodeVec);
	for (size_t i = 0; i < nodeVec.size(); i++)
	{
		size_t index = nodeVec[i].getindex();
		TriNodeMap_[index] = nodeVec[i];
		cgalPointsIndexVec_.emplace_back(index);
	}
}
void TriBase::insertCtrlLine(std::map<size_t, std::vector<double>>CoorMap)
{
	for (const auto& it : CoorMap)
	{
		size_t index = it.first;
		std::vector<double> coor = it.second;
		double x = 0, y = 0, z = 0;
		if (coor.size() == 2) {
			x = coor[0], y = coor[1];
		}
		else if (coor.size() == 3) {
			x = coor[0], y = coor[1], z = coor[2];
		}
		TriNodeMap_[index] = TriNode(x, y, z, index, TriNodeType::inner);
		cgalPointsIndexVec_.emplace_back(index);
	}
}
void TriBase::insertPoint(std::map<size_t, std::vector<double>>CoorMap)
{
	for (const auto& it : CoorMap)
	{
		size_t index = it.first;
		std::vector<double> coor = it.second;
		double x = 0, y = 0, z = 0;
		if (coor.size() == 2) {
			x = coor[0], y = coor[1];
		}
		else if (coor.size() == 3) {
			x = coor[0], y = coor[1], z = coor[2];
		}
		TriNodeMap_[index] = TriNode(x, y, z, index);
		cgalPointsIndexVec_.emplace_back(index);
	}
}


/*----------------------------------------------------------------------------------------------
	本部分是一个简单的内凹外凹多边形；用于DEBUG
------------------------------------------------------------------------------------------------*/
TriBase::TriBase(const std::string& type)
{
	std::vector<std::vector<double>> poly_1_outter = { {1,2}, {3,2}, {7,2}, {7,4}, {8,7}, {5,8}, {2,6}, {3,4} };
	std::vector<std::vector<double>> poly_2_inner = { {4,4}, {6,3}, {6,6}, {4,6}, {5,5} };
	size_t index = 1;
	for (const auto& per_it : poly_1_outter)
	{
		TriNode Point_cdt{ per_it[0],per_it[1],0, index, TriNodeType::outter };
		TriNodeMap_[index] = Point_cdt;
		cgalPointsIndexVec_.emplace_back(index);
		index++;
	}
	for (const auto& per_it : poly_2_inner)
	{
		TriNode Point_cdt{ per_it[0],per_it[1],0, index, TriNodeType::inner };
		TriNodeMap_[index] = Point_cdt;
		cgalPointsIndexVec_.emplace_back(index);
		index++;
	}
	dividType_ = type;
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
	std::string basePath = "tempData";
	createFolder(basePath);
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
		tecplotFile << p.x() << " " << p.y() << "\n";
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
	std::vector<TriNode> poly_outerCoor;
	std::vector<TriNode> poly_innerCoor;
	//TriNodeMap_是没有顺序的，在TriNodeMap_中循环就会导致控制边界也没有顺序了，插入边界时候出错。
	//在cgalPointsIndexVec_中循环就保证了顺序。
	for (const auto& perit : cgalPointsIndexVec_)
	{
		const auto& cur_trinode = TriNodeMap_.at(perit);
		if (cur_trinode.gettype() == TriNodeType::outter)
		{
			cgal_outerPointsVec.push_back({ cur_trinode.x() ,cur_trinode.y() });
			poly_outerCoor.emplace_back(cur_trinode);
		}
		else if (cur_trinode.gettype() == TriNodeType::inner)
		{
			cgal_innerPointsVec.push_back({ cur_trinode.x() ,cur_trinode.y() });
			poly_innerCoor.emplace_back(cur_trinode);
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
		for (size_t i = 0; i < cgal_innerPointsVec.size() - 1; ++i)
			cdt.insert_constraint(cgal_innerPointsVec[i], cgal_innerPointsVec[i + 1]);
	}

	//step 6: 写入三角形索引
	size_t triIndex = 0;
	for (Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit)
	{
		size_t v0 = vertex_to_index[fit->vertex(0)];
		size_t v1 = vertex_to_index[fit->vertex(1)];
		size_t v2 = vertex_to_index[fit->vertex(2)];
		TriNode p0 = index2Coor(v0);
		TriNode p1 = index2Coor(v1);
		TriNode p2 = index2Coor(v2);
		//step 1 : 判断三角形的三个顶点是否在一条直线上,
		if (isCollinear(p0, p1, p2))
			continue;
		TriNode triCenter = (p0 += p1 += p2) /= 3.0;
		if (dividType_ == "ring_ring") {
			//step 2 : 判断三角形是否在<内层>边界多边形<内部>，如果<在内部>则删除
			if (isPointInPolygon_moreInner(triCenter, poly_innerCoor))
				continue;
			//step 3 : 判断三角形是否在<外层>边界多边形<外部>，如果<在外部>则删除
			if (!isPointInPolygon_moreOuter(triCenter, poly_outerCoor))
				continue;
		}
		else if (dividType_ == "ring_line") {
			if (!isPointInPolygon_moreOuter(triCenter, poly_outerCoor))
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
			for (auto& var : triangles)
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
	std::vector<std::pair<TriNode, size_t >> newTri_1_pair;
	std::vector<std::pair<TriNode, size_t >> newTri_2_pair;
	TriNode average_1;
	TriNode average_2;
	for (auto it : newTri_1)
	{
		TriNode temcoor = index2Coor(it);
		newTri_1_pair.push_back(std::make_pair(temcoor, it));
		average_1 += (temcoor /= (double)newTri_1.size());
	}
	for (auto it : newTri_2)
	{
		TriNode temcoor = index2Coor(it);
		newTri_2_pair.push_back(std::make_pair(temcoor, it));
		average_2 += (temcoor /= (double)newTri_2.size());
	}
	reorderPointsIndex(newTri_1_pair, average_1);//按逆时针重新排序
	reorderPointsIndex(newTri_2_pair, average_2);
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
	addEdgeTableItem(current_triptr);

	std::vector<TriNode> triNodes_2;
	triNodes_2.emplace_back(TriNodeMap_.at(newTri_2.at(0)));
	triNodes_2.emplace_back(TriNodeMap_.at(newTri_2.at(1)));
	triNodes_2.emplace_back(TriNodeMap_.at(newTri_2.at(2)));
	TriEle triEle_2(triNodes_2, GetTriEleVec().size());
	triEle_2.setTriCenter(average_2);
	const auto& another_ptr = std::make_shared<TriEle>(triEle_2);
	addTriEle(another_ptr);
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
	std::vector<std::pair<TriNode, size_t >> quad_pair;
	TriNode average;
	for (auto it : quad)
	{
		TriNode temcoor = index2Coor(it);
		quad_pair.push_back(std::make_pair(temcoor, it));
		average += (temcoor /= (double)quad.size());
	}

	reorderPointsIndex(quad_pair, average);//按逆时针重新排序
	quad.clear();
	for (auto& it : quad_pair)
		quad.push_back(it.second);

	getCell().push_back(quad);//第一次getCell是cell_是空的，在此添加cell_

	//标记相邻三角形为TriEleTag::used
	for (auto& var : triangles)
	{
		if (var->getTriEleTag() == TriEleTag::used) {
			//spdlog::error("索引为{}的三角形转换四边形出错!", var->GetTriEleIndex());
			throw std::runtime_error("Caught an exception!");
		}
		var->setTriEleTag(TriEleTag::used);
	}

	return;
}

void TriBase::gencell(const triEleMap_& content)
{
	for (auto& triEle : content)
		cell_.push_back(triEle.second->GetTriEleVertex());
	return;
}

void TriBase::boostMesh()
{
	genCDT();
	//reconstructNotIdealTriGrid();//:Q重构非理想三角形是必须的吗？A:不是
	std::string method = "CDT";
	//method = "HYB";//用于DEBUG
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

//方法函数--------------------------------------------------------------------------------------------------------------------------
void TriBase::reorderPointsIndex(std::vector<std::pair<TriNode, size_t>>& PointLists, const TriNode& modelCenter)
{
	std::size_t wallNodeNum = PointLists.size();
	std::vector<std::pair<double, size_t>> angleAndIndex;
	angleAndIndex.resize(wallNodeNum);
	for (size_t iNode = 0; iNode < wallNodeNum; ++iNode)
	{
		auto& perNode = PointLists.at(iNode).first;
		TriNode vec = perNode - modelCenter;
		double angle = atan2(vec.y(), vec.x());
		angleAndIndex[iNode] = std::make_pair(angle, iNode);
	}
	std::sort(angleAndIndex.begin(), angleAndIndex.end(),
		[](const std::pair<double, size_t>& a, const std::pair<double, size_t>& b) {return a.first < b.first; });

	std::vector<size_t> reorderedIndices;
	std::vector<TriNode> NodeVec_temp;
	std::vector<std::pair<TriNode, size_t>> NodePairVec_temp2;
	NodeVec_temp.resize(wallNodeNum);
	reorderedIndices.resize(wallNodeNum);
	for (size_t i = 0; i < wallNodeNum; ++i)
	{
		NodeVec_temp[i] = PointLists[angleAndIndex[i].second].first;
		reorderedIndices[i] = PointLists[angleAndIndex[i].second].second;
		NodePairVec_temp2.push_back(std::make_pair(NodeVec_temp[i], reorderedIndices[i]));
	}
	PointLists = NodePairVec_temp2;//会修改原本点的顺序;
}

bool TriBase::isCollinear(const TriNode& p1, const TriNode& p2, const TriNode& p3)
{
	// 使用斜率法判断三点是否共线  
	// 如果斜率相同或任意两点重合（在容差范围内），则它们共线  
	constexpr double epsilon = std::numeric_limits<double>::epsilon();
	double dx1 = p2.x() - p1.x();
	double dy1 = p2.y() - p1.y();
	double dx2 = p3.x() - p1.x();
	double dy2 = p3.y() - p1.y();

	// 使用斜率来判断是否共线  
	if (std::abs(dx1 * dy2 - dx2 * dy1) < epsilon) {
		return true;
	}
	// 检查是否接近0（即检查是否有重合的点）  
	if (std::abs(dx1) < epsilon && std::abs(dx2) < epsilon) {
		// 如果x坐标都接近相同，则检查y坐标是否也接近相同  
		return std::abs(dy1) < epsilon && std::abs(dy2) < epsilon;
	}

	return false;
}

void TriBase::createFolder(const std::string& basePath)const
{
	// 构造完整的文件夹路径  
	std::filesystem::path folderPath = basePath;

	// 检查文件夹是否已经存在  
	if (!std::filesystem::exists(folderPath)) {
		// 创建文件夹  
		std::filesystem::create_directories(folderPath);
		//ZaranLog::info("Folder '{}' created successfully!", basePath);
	}
}

/*-------------------------------------------------------------------
*  Function: isPointInPolygon_moreInner
*  Purpose: 判断点是否在多边形内部版本1：包括与边界重合或顶点重合的点都认为在内部
*  Arguments:
*    p - 目标判断点
*    polygon - 有序封闭多边形
*  Returns:
*    boll
-------------------------------------------------------------------*/
bool TriBase::isPointInPolygon_moreInner(const TriNode& p, const std::vector<TriNode>& polygon) {
	const double EPSILON = 1e-10;
	bool inside = false;
	size_t j = polygon.size() - 1;

	for (size_t i = 0; i < polygon.size(); i++) {
		// Check if point is on vertex
		if ((std::fabs(p.x() - polygon[i].x()) < EPSILON) && (std::fabs(p.y() - polygon[i].y()) < EPSILON)) {
			return true; // Point is on a vertex
		}

		// Check if point is on edge
		double minX = std::min(polygon[i].x(), polygon[j].x());
		double maxX = std::max(polygon[i].x(), polygon[j].x());
		double minY = std::min(polygon[i].y(), polygon[j].y());
		double maxY = std::max(polygon[i].y(), polygon[j].y());

		if ((minX <= p.x() && p.x() <= maxX) && (minY <= p.y() && p.y() <= maxY)) {
			double dx = polygon[j].x() - polygon[i].x();
			double dy = polygon[j].y() - polygon[i].y();
			if (std::fabs(dy) < EPSILON) { // horizontal line
				if (std::fabs(p.y() - polygon[i].y()) < EPSILON) {
					return true; // Point is on horizontal edge
				}
			}
			else if (std::fabs(dx) < EPSILON) { // vertical line
				if (std::fabs(p.x() - polygon[i].x()) < EPSILON) {
					return true; // Point is on vertical edge
				}
			}
			else {
				double slope = dy / dx;
				double intercept = polygon[i].y() - slope * polygon[i].x();
				if (std::fabs(p.y() - (slope * p.x() + intercept)) < EPSILON) {
					return true; // Point is on non-vertical edge
				}
			}
		}

		// Check if point is inside using ray-casting algorithm
		if (((polygon[i].y() <= p.y()) && (p.y() < polygon[j].y())) ||
			((polygon[j].y() <= p.y()) && (p.y() < polygon[i].y()))) {
			double intersectionX = polygon[i].x() + (p.y() - polygon[i].y()) * (polygon[j].x() - polygon[i].x()) / (polygon[j].y() - polygon[i].y());
			if (intersectionX < p.x()) {
				inside = !inside;
			}
		}

		j = i;
	}
	return inside;
}

/*-------------------------------------------------------------------
*  Function: isPointInPolygon_moreOuter
*  Purpose: 判断点是否在多边形内部版本2：包括与边界重合或顶点重合的点都认为在外部
*  Arguments:
*    p - 目标判断点
*    polygon - 有序封闭多边形
*  Returns:
*    boll
-------------------------------------------------------------------*/
bool TriBase::isPointInPolygon_moreOuter(const TriNode& p, const std::vector<TriNode>& polygon) {
	const double EPSILON = 1e-10;
	bool inside = false;
	size_t j = polygon.size() - 1;

	for (size_t i = 0; i < polygon.size(); i++) {
		// Check if point is on vertex
		if ((std::fabs(p.x() - polygon[i].x()) < EPSILON) && (std::fabs(p.y() - polygon[i].y()) < EPSILON)) {
			return false; // Point is on a vertex
		}

		// Check if point is on edge
		double minX = std::min(polygon[i].x(), polygon[j].x());
		double maxX = std::max(polygon[i].x(), polygon[j].x());
		double minY = std::min(polygon[i].y(), polygon[j].y());
		double maxY = std::max(polygon[i].y(), polygon[j].y());

		if ((minX <= p.x() && p.x() <= maxX) && (minY <= p.y() && p.y() <= maxY)) {
			double dx = polygon[j].x() - polygon[i].x();
			double dy = polygon[j].y() - polygon[i].y();
			if (std::fabs(dy) < EPSILON) { // horizontal line
				if (std::fabs(p.y() - polygon[i].y()) < EPSILON) {
					return false; // Point is on horizontal edge
				}
			}
			else if (std::fabs(dx) < EPSILON) { // vertical line
				if (std::fabs(p.x() - polygon[i].x()) < EPSILON) {
					return false; // Point is on vertical edge
				}
			}
			else {
				double slope = dy / dx;
				double intercept = polygon[i].y() - slope * polygon[i].x();
				if (std::fabs(p.y() - (slope * p.x() + intercept)) < EPSILON) {
					return false; // Point is on non-vertical edge
				}
			}
		}

		// Check if point is inside using ray-casting algorithm
		if (((polygon[i].y() <= p.y()) && (p.y() < polygon[j].y())) ||
			((polygon[j].y() <= p.y()) && (p.y() < polygon[i].y()))) {
			double intersectionX = polygon[i].x() + (p.y() - polygon[i].y()) * (polygon[j].x() - polygon[i].x()) / (polygon[j].y() - polygon[i].y());
			if (intersectionX < p.x()) {
				inside = !inside;
			}
		}

		j = i;
	}
	return inside;
}

void TriBase::reorderNodes(std::vector<TriNode>& elements)
{
	//step 1: 计算中心点
	std::size_t elementNum = elements.size();
	TriNode modelCenter{ 0,0,0 };
	for (size_t iPoint = 0; iPoint < elementNum; ++iPoint)
	{
		TriNode per = elements[iPoint];
		modelCenter += per;
	}
	modelCenter /= (double)elementNum;

	//step 2: 计算角度并排序
	std::vector<std::pair<double, size_t>> angleAndIndex;
	angleAndIndex.reserve(elementNum);
	for (size_t i = 0; i < elementNum; ++i)
	{
		auto& element = elements[i];
		TriNode vec = element - modelCenter;
		double angle = atan2(vec.y(), vec.x());
		angleAndIndex.emplace_back(angle, i);
	}

	std::sort(angleAndIndex.begin(), angleAndIndex.end(),
		[](const std::pair<double, size_t>& a, const std::pair<double, size_t>& b) { return a.first < b.first; });

	std::vector<TriNode> elementsTemp(elementNum);
	for (size_t i = 0; i < elementNum; ++i)
	{
		elementsTemp[i] = elements[angleAndIndex[i].second];
	}

	elements = std::move(elementsTemp); // 通过移动语义优化性能
}