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
		if (edge->vA == edge->vB) {
			//spdlog::error("��ǰ�ߵ�TriNode����vA={},vB={}", edge->vA, edge->vB);
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

	//����TriNodeMap��ÿ��TriNode������ʵ������ʱ��͸���������,ÿ�������������Ψһ��
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
	��������һ���򵥵��ڰ��ⰼ����Σ�����DEBUG
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
	std::string basePath = "tempData";
	createFolder(basePath);
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
		tecplotFile << p.x() << " " << p.y() << "\n";
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
	std::vector<TriNode> poly_outerCoor;
	std::vector<TriNode> poly_innerCoor;
	//TriNodeMap_��û��˳��ģ���TriNodeMap_��ѭ���ͻᵼ�¿��Ʊ߽�Ҳû��˳���ˣ�����߽�ʱ�����
	//��cgalPointsIndexVec_��ѭ���ͱ�֤��˳��
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
		for (size_t i = 0; i < cgal_innerPointsVec.size() - 1; ++i)
			cdt.insert_constraint(cgal_innerPointsVec[i], cgal_innerPointsVec[i + 1]);
	}

	//step 6: д������������
	size_t triIndex = 0;
	for (Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit)
	{
		size_t v0 = vertex_to_index[fit->vertex(0)];
		size_t v1 = vertex_to_index[fit->vertex(1)];
		size_t v2 = vertex_to_index[fit->vertex(2)];
		TriNode p0 = index2Coor(v0);
		TriNode p1 = index2Coor(v1);
		TriNode p2 = index2Coor(v2);
		//step 1 : �ж������ε����������Ƿ���һ��ֱ����,
		if (isCollinear(p0, p1, p2))
			continue;
		TriNode triCenter = (p0 += p1 += p2) /= 3.0;
		if (dividType_ == "ring_ring") {
			//step 2 : �ж��������Ƿ���<�ڲ�>�߽�����<�ڲ�>�����<���ڲ�>��ɾ��
			if (isPointInPolygon_moreInner(triCenter, poly_innerCoor))
				continue;
			//step 3 : �ж��������Ƿ���<���>�߽�����<�ⲿ>�����<���ⲿ>��ɾ��
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
			for (auto& var : triangles)
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
	reorderPointsIndex(newTri_1_pair, average_1);//����ʱ����������
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
	std::vector<std::pair<TriNode, size_t >> quad_pair;
	TriNode average;
	for (auto it : quad)
	{
		TriNode temcoor = index2Coor(it);
		quad_pair.push_back(std::make_pair(temcoor, it));
		average += (temcoor /= (double)quad.size());
	}

	reorderPointsIndex(quad_pair, average);//����ʱ����������
	quad.clear();
	for (auto& it : quad_pair)
		quad.push_back(it.second);

	getCell().push_back(quad);//��һ��getCell��cell_�ǿյģ��ڴ����cell_

	//�������������ΪTriEleTag::used
	for (auto& var : triangles)
	{
		if (var->getTriEleTag() == TriEleTag::used) {
			//spdlog::error("����Ϊ{}��������ת���ı��γ���!", var->GetTriEleIndex());
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
	//reconstructNotIdealTriGrid();//:Q�ع��������������Ǳ������A:����
	std::string method = "CDT";
	//method = "HYB";//����DEBUG
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

//��������--------------------------------------------------------------------------------------------------------------------------
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
	PointLists = NodePairVec_temp2;//���޸�ԭ�����˳��;
}

bool TriBase::isCollinear(const TriNode& p1, const TriNode& p2, const TriNode& p3)
{
	// ʹ��б�ʷ��ж������Ƿ���  
	// ���б����ͬ�����������غϣ����ݲΧ�ڣ��������ǹ���  
	constexpr double epsilon = std::numeric_limits<double>::epsilon();
	double dx1 = p2.x() - p1.x();
	double dy1 = p2.y() - p1.y();
	double dx2 = p3.x() - p1.x();
	double dy2 = p3.y() - p1.y();

	// ʹ��б�����ж��Ƿ���  
	if (std::abs(dx1 * dy2 - dx2 * dy1) < epsilon) {
		return true;
	}
	// ����Ƿ�ӽ�0��������Ƿ����غϵĵ㣩  
	if (std::abs(dx1) < epsilon && std::abs(dx2) < epsilon) {
		// ���x���궼�ӽ���ͬ������y�����Ƿ�Ҳ�ӽ���ͬ  
		return std::abs(dy1) < epsilon && std::abs(dy2) < epsilon;
	}

	return false;
}

void TriBase::createFolder(const std::string& basePath)const
{
	// �����������ļ���·��  
	std::filesystem::path folderPath = basePath;

	// ����ļ����Ƿ��Ѿ�����  
	if (!std::filesystem::exists(folderPath)) {
		// �����ļ���  
		std::filesystem::create_directories(folderPath);
		//ZaranLog::info("Folder '{}' created successfully!", basePath);
	}
}

/*-------------------------------------------------------------------
*  Function: isPointInPolygon_moreInner
*  Purpose: �жϵ��Ƿ��ڶ�����ڲ��汾1��������߽��غϻ򶥵��غϵĵ㶼��Ϊ���ڲ�
*  Arguments:
*    p - Ŀ���жϵ�
*    polygon - �����ն����
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
*  Purpose: �жϵ��Ƿ��ڶ�����ڲ��汾2��������߽��غϻ򶥵��غϵĵ㶼��Ϊ���ⲿ
*  Arguments:
*    p - Ŀ���жϵ�
*    polygon - �����ն����
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
	//step 1: �������ĵ�
	std::size_t elementNum = elements.size();
	TriNode modelCenter{ 0,0,0 };
	for (size_t iPoint = 0; iPoint < elementNum; ++iPoint)
	{
		TriNode per = elements[iPoint];
		modelCenter += per;
	}
	modelCenter /= (double)elementNum;

	//step 2: ����ǶȲ�����
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

	elements = std::move(elementsTemp); // ͨ���ƶ������Ż�����
}