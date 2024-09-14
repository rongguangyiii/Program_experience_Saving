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
*      ��������ļ�(TriEdge),ͨ��FNV-1a��ϣ�㷨����õ�һ����ϣֵ(Ψһ����)
*  Description: 
*     std::unordered_map ���Զ�ʹ��EdgeHash��operator()������edge�Ĺ�ϣֵ�����ݴ�
*     ȷ���ڹ�ϣ����Ӧ�ô洢�����ֵ�Ե�λ�á��㲻��Ҫ��Ҳ��Ӧ�ã�ֱ�ӵ���EdgeHash��operator()��
*     ����std::unordered_map���ڲ�ʵ��ϸ�ڣ�������ȷ������Ҫʱ��ȷ��ʹ�ù�ϣ������
*  Usage: 
*		std::unordered_map<TriEdge, std::vector<TriEle>, EdgeHash>
*  Notice : 
*		��Ҫ���� ��(TriEdge) ��Ӧ��operator==()��������hashֵʵ�ֺ�����Ҳ���Ǳ�������
*  Arguments:
*    k - ����ļ������ߣ��������������
*  Returns:
*   size_t - ����һ����ϣֵ
---------------------------------------------------------------------------------------------*/
std::size_t EdgeHash::operator()(const std::shared_ptr<TriEdge>& k) const {
	// ʹ��FNV-1a��ϣ�㷨�ı��壬��������������ͨ����������  
	std::size_t hash = 2166136261; // FNV offset basis  
	hash ^= k->vA;
	hash *= 16777619; // FNV prime  
	hash ^= k->vB;
	hash *= 16777619; // �ٴλ��  
	return hash;
}

void TriBase::addEdgeTableItem(const triElePtr_& triEle_ptr)
{
	// ��ÿ���ߣ�������ӵ��߱���
	for (const auto& edge : triEle_ptr->getTriEdgesRef()) 
	{
		if (edge->vA == edge->vB){
			//spdlog::error("��ǰ�ߵ�TriNode����vA={},vB={}", edge->vA, edge->vB);
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

	//����TriNodeMap��ÿ��TriNode������ʵ������ʱ��͸���������,ÿ�������������Ψһ��
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
	��������һ���򵥵��ڰ��ⰼ����Σ�����DEBUG
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
*  Purpose: ���±߱�(��ϣֵ��)edgeTable_����һ��û�б߱�Ļ��ʹ�����ϣֵ��
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
*  Purpose: ����cell_,���ڴ洢�����κ��ı��εĶ�������
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
*  Purpose: ��������κ��ı������� Tecplot �ļ��У����ڲ鿴�����
*  Arguments:
*    pointIndexVec - ���еĵ��������
*    cellVec - ���е������κ�˼�棬�������Tecplot��ԪԪ�����ӱ�
*    vertex2tecplotIndex - ����ӳ��map,����ת���ڵ�������0��ʼ��
*    title - Tecplot�ļ�����
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
	// ����һ�� Tecplot �ļ�

	std::ofstream tecplotFile(basePath + "/" + outfilename);
	tecplotFile << "TITLE = \"Triangle Output   LGYYYY\"\n";
	tecplotFile << "VARIABLES = \"X\", \"Y\"\n";

	// д�붥������
	tecplotFile << "ZONE T=\"Triangles\", N=" << pointIndexVec.size()
		<< ", E=" << cellVec.size() << ", F=FEPOINT, ET=QUADRILATERAL\n";
	for (const auto& vit : pointIndexVec)
	{
		const auto& p = index2Coor(vit);
		tecplotFile << p.x << " " << p.y << "\n";
		//spdlog::debug("����ڵ����꣺{},{}", p.x(), p.y());
	}

	// д�뵥Ԫ��Ϣ
	for (const auto& cell : cellVec)
	{
		if (cell.size() == 3) { // �����ε�Ԫ
			tecplotFile << vertex2tecplotIndex.at(cell[0]) + 1 << " " << vertex2tecplotIndex.at(cell[1]) + 1 << " "
				<< vertex2tecplotIndex.at(cell[2]) + 1 << " " << vertex2tecplotIndex.at(cell[0]) + 1 << std::endl;
		}
		else if (cell.size() == 4) { // �ı��ε�Ԫ
			tecplotFile << vertex2tecplotIndex.at(cell[0]) + 1 << " " << vertex2tecplotIndex.at(cell[1]) + 1 << " "
				<< vertex2tecplotIndex.at(cell[2]) + 1 << " " << vertex2tecplotIndex.at(cell[3]) + 1 << std::endl;
		}
		else {
			//spdlog::error("������OutputConnectGrid2Tecplot(): ��Ԫ�Ķ����������ԣ�����");
		}
	}

	tecplotFile.close();

	return;
}

/*-------------------------------------------------------------------
*  Function: Gen_initialCDT
*  Purpose: ���ɳ�ʼ���������������Tecplot�ļ��У����ڲ鿴�����
*  Arguments:
*    null
*  Returns:
*    null
-------------------------------------------------------------------*/
void TriBase::genCDT()
{
	//step 1: ����CDTԼ���߽��������ɸѡ�߽硣
	std::vector<Point_2> cgal_outerPointsVec;
	std::vector<Point_2> cgal_innerPointsVec;
	std::vector<Coord> poly_outerCoor;
	std::vector<Coord> poly_innerCoor;
	//TriNodeMap_��û��˳��ģ���TriNodeMap_��ѭ���ͻᵼ�¿��Ʊ߽�Ҳû��˳���ˣ�����߽�ʱ�����
	//��cgalPointsIndexVec_��ѭ���ͱ�֤��˳��
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

	//step 2: �ϲ������㼯  cgalPointsIndexVec_��˳�������cgal_PointsVec��˳�򣡣���
	std::vector<Point_2> cgal_PointsVec(cgal_outerPointsVec.begin(), cgal_outerPointsVec.end());
	cgal_PointsVec.insert(cgal_PointsVec.end(), cgal_innerPointsVec.begin(), cgal_innerPointsVec.end());

	//step 3: ����һ��Լ�������ʷ�
	CDT cdt;

	//step 4: �����������������
	std::unordered_map<size_t, size_t> vertex2tecplotIndex;//���ڵ�����ӳ�䵽�ֲ���0��ʼ���������������Tecplot�ļ�
	std::map<Vertex_handle, size_t> vertex_to_index;//��CDT����ӳ�䵽�ڵ�����
	for (size_t iPoint = 0; iPoint != cgalPointsIndexVec_.size(); ++iPoint)
	{
		Vertex_handle vh = cdt.insert(cgal_PointsVec.at(iPoint));
		vertex_to_index[vh] = cgalPointsIndexVec_.at(iPoint);
		vertex2tecplotIndex[cgalPointsIndexVec_.at(iPoint)] = iPoint;
	}
	vertex2tecplotIndex_.first = cgalPointsIndexVec_;//�ڵ�����
	vertex2tecplotIndex_.second = vertex2tecplotIndex;//�ֲ�����
	
	/*-----------------------------------------------------------------
	 * Note��CGAL���insert_constraint���������ڲ���Լ���߽磻
	 *	     ͨ����Ҫ��������ΪPoint_2�Ĳ�������ʾ����Լ���߶ε������˵㣬
	 *       ���ж���߶β���ʱ��ѡ����λ�����Ľڵ����������ġ�
	 *-----------------------------------------------------------------*/
	//step 5: ����Լ���߽�,��Ϊ������Բ��Լ�������߽�Ҫ�ֱ���룡����
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

	//step 6: д������������
	size_t triIndex = 0;
	for (Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit)
	{
		size_t v0 = vertex_to_index[fit->vertex(0)];
		size_t v1 = vertex_to_index[fit->vertex(1)];
		size_t v2 = vertex_to_index[fit->vertex(2)];
		Coord p0 = index2Coor(v0);
		Coord p1 = index2Coor(v1);
		Coord p2 = index2Coor(v2);
		//step 1 : �ж������ε����������Ƿ���һ��ֱ����,
		if (Tools::isCollinear(p0, p1, p2))
			continue;
		Coord triCenter = (p0 += p1 += p2) /= 3.0;
		if (dividType_ == "ring_ring")
		{
			//step 3 : �ж��������Ƿ���<���>�߽�����<�ⲿ>�����<���ⲿ>��ɾ��
			if (!Tools::isPointInPolygon_moreOuter(triCenter, poly_outerCoor))
				continue;
			//step 2 : �ж��������Ƿ���<�ڲ�>�߽�����<�ڲ�>�����<���ڲ�>��ɾ��
			if (Tools::isPointInPolygon_moreInner(triCenter, poly_innerCoor))
				continue;
		}
		else if (dividType_ == "ring_points") {
			//step 2 : �ж��������Ƿ���<�ڲ�>�߽�����<�ڲ�>�����<���ڲ�>��ɾ��
			if (Tools::isPointInPolygon_moreInner(triCenter, poly_innerCoor))
				continue;
		}
		else if (dividType_ == "ring_line") {
			//step 3 : �ж��������Ƿ���<���>�߽�����<�ⲿ>�����<���ⲿ>��ɾ��
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

	//step 7: �����ʼ����������
	std::vector< std::vector<size_t> > cellVec;
	for (const auto& perTri : triEleVec_)
		cellVec.push_back(perTri.second->GetTriEleVertex());
	OutputConnectGrid2Tecplot(vertex2tecplotIndex_.first, cellVec, vertex2tecplotIndex_.second, "Tri_ini");

	//step ... ���������Դ�����������Ϣ��

	return;
}

/*-------------------------------------------------------------------
*  Function: filterTriGrid
*  Purpose: ɸѡ����������ɾ�����ڲ���������ڵĺ��ⲿ���������������Ρ�
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
		//3�����㶼���ⲿ�߽��ϻ��ڲ��߽���
		if (currentTri->getTriEleType() == TriEleType::Inner_only ||
			currentTri->getTriEleType() == TriEleType::Outer_only) {
			addNotIdealTri(*currentTri);
		}
	}
	return;
}

/*-------------------------------------------------------------------
*  Function: reconstructNotIdealTriGrid
*  Purpose: �Բ�����Ҫ��ģ������ڻ����⻷֮��������ν��д���
*			����������������������
*  Arguments:
*    null
*  Returns:
*    null
-------------------------------------------------------------------*/
void TriBase::reconstructNotIdealTriGrid()
{
	//step 0: ɸѡ��Ҫ�ع�������������
	filterTri();
	std::vector<TriEle>& notIdealTriVec = getNotIdealTriVec();
	//step 0:���û�в�����������Σ�ֱ�ӷ���
	if (0 == notIdealTriVec.size())
		return;

	//step 1: ���ɱ߱�
	genEdgeTable(); 
	//step 2: �ع�������������
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
			//Ϊ��ֹ��ɾ�����Σ��Ƚ�Ҫɾ���������������洢����
			for(auto& var : triangles)
				toDeleteTriIndex.push_back(var->GetTriEleIndex());
			reconstructNotIdealTriGrid(*edge, triangles);
		}
	}
	//step 3: ɾ�����������β�ͬʱ���±߱�
	for (size_t cc : toDeleteTriIndex)
		removeTriEle(cc);

	//step 4: �������
	std::vector< std::vector<size_t> > cellVec;
	for (const auto& perTri : GetTriEleVec())
		cellVec.push_back(perTri.second->GetTriEleVertex());
	OutputConnectGrid2Tecplot(vertex2tecplotIndex_.first, cellVec, vertex2tecplotIndex_.second, "Tri_reconstruct");
	return;
}

/*-------------------------------------------------------------------
*  Function: reconstructNotIdealTriGrid(edge,triangles)
*  Purpose: �Ե�ǰ���ߵ������������������ӡ�
*  Arguments:
*    edge - ����
*   triangles - ���ߵ�����������
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
	//�ҵ���������coindex_a��coindex_b��ͬ�Ķ���
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
	Tools::reorderPointsIndex(newTri_1_pair, average_1);//����ʱ����������
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
*  Purpose: ���ɱ�������ͷǽṹ��������Ӳ��֡����ò���CDT������ı����������
*  Arguments:
*    null
*  Returns:
*    null
-------------------------------------------------------------------*/
void TriBase::genTriAndQuad()
{
	//step 0 :�ж��Ƿ�Ҫ���ɱ߱�
	if (0 == edgeTable_.size())
		genEdgeTable();

	//step 1 :������תΪ�ı���
	for (const auto& perit : edgeTable_)
	{
		auto& cur_edge = perit.first;
		auto& cur_triVec = perit.second;
		if (cur_triVec.size() == 0 || cur_triVec.size() >= 3) {
			//spdlog::error("�߱������������GenTriAndQuad() 835");
			throw std::runtime_error("Caught an exception!");
		}
		//��ǰ���Ǳ߽�ߡ�
		if (cur_triVec.size() == 1)
			continue;
		//��ǰ�߶�Ӧ������������������һ���Ѿ��������ӹ���
		if (cur_triVec[0]->getTriEleTag() == TriEleTag::used ||
			cur_triVec[1]->getTriEleTag() == TriEleTag::used)
			continue;
		//��ǰ�߶�Ӧ���������������Ͷ���һ���ġ�
		if (cur_triVec[0]->getTriEleType() == cur_triVec[1]->getTriEleType())
			continue;

		//�����ǰ�������з���Ҫ������������Σ��ͽ����ع�Ϊ�ı��Σ�����ͬʱ�������������ΪTriEleTag::used��
		genTriAndQuad(*cur_edge, cur_triVec);

	}//end for edge

	//step 2 :����������ѭ���������ٴ�����Ȼδ��ʹ�õ�������
	GenCompleteCell();
	//step 3 :�������
	OutputConnectGrid2Tecplot(vertex2tecplotIndex_.first, getCell(), vertex2tecplotIndex_.second, "Tri_quad_final");

	return;
}

/*-------------------------------------------------------------------
*  Function: GenUnstructConnectElement_Tri2quad
*  Purpose: ���ɱ�������ͷǽṹ��������Ӳ��֡����ò���CDT������ı����������
*  Arguments:
*    edge - ����
*   triangles - ���ߵ�����������
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
	//�ҵ���������coindex_a��coindex_b��ͬ�Ķ���
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

	Tools::reorderPointsIndex(quad_pair, average);//����ʱ����������
	quad.clear();
	for (auto& it : quad_pair)
		quad.push_back(it.second);

	getCell().push_back(quad);//��һ��getCell��cell_�ǿյģ��ڴ����cell_

	//�������������ΪTriEleTag::used
	for (auto& var : triangles)
	{
		if (var->getTriEleTag() == TriEleTag::used){
			//spdlog::error("����Ϊ{}��������ת���ı��γ���!", var->GetTriEleIndex());
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
	//reconstructNotIdealTriGrid();//:Q�ع��������������Ǳ������A�о����Ǵ�������֤...
	//std::string method = GlobalData::GetString("ConnectMethod");
	std::string method = "HYB";//����DEBUG
	if (method == "CDT")
		gencell(GetTriEleVec());
	else if (method == "HYB")
		genTriAndQuad();
	else
		//spdlog::error("gencell()����key word error��");
	return;
}

void TriBase::addTriEle(const triElePtr_& triEle)
{ 
	size_t index = triEle->GetTriEleIndex();
	triEleVec_.insert_or_assign(index, triEle); // ����Ԫ���Ƿ��Ѿ����ڣ���������ȷ�������¡�
	//triEleVec_.insert_or_assign(index, std::move(triEle)); // ����Ԫ���Ƿ��Ѿ����ڣ���������ȷ�������¡�
};

void TriBase::removeTriFromEdgeTable(const triElePtr_& triEle_ptr) {
	// ��ÿ���ߣ��ӱ߱����Ƴ�triEle_ptr
	for (const auto& edge : triEle_ptr->getTriEdgesRef()) {
		auto it = edgeTable_.find(edge);
		if (it != edgeTable_.end()) {
			auto& triEleList = it->second;
			triEleList.erase(std::remove(triEleList.begin(), triEleList.end(), triEle_ptr), triEleList.end());
			if (triEleList.empty()) {
				edgeTable_.erase(it); // ����б�Ϊ�գ��Ƴ��ñ���Ŀ
			}
		}
	}
}

void TriBase::removeTriEle(size_t triEleKey) {
	auto it = triEleVec_.find(triEleKey);
	if (it != triEleVec_.end()) {
		auto& triEle_ptr = it->second; // ʹ�����еĹ���ָ��
		removeTriFromEdgeTable(triEle_ptr); // ���ȴӱ߱����Ƴ�
		triEleVec_.erase(it); // Ȼ��� triEleVec_ ���Ƴ�
	}
	else {
		throw std::runtime_error("TriEle not found in triEleVec_");
	}
}