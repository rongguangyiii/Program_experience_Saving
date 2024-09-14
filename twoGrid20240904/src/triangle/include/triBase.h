/*---------------------------------------------------------------------------
	ZaRan	-	A Totally Automatic CFD Software
	Copyright (C) ,Since 2020
-------------------------------------------------------------------------------
License
	This file is part of ZaRan.

!	@file	tri_info.h
!	@brief	Generate ZaRan Grid.
!	@author	Liu Guangying.
!   @date  2024.04.20
\*---------------------------------------------------------------------------*/
#ifndef TRI_ELEBASE_H  
#define TRI_ELEBASE_H  
#pragma once
#include "triangle/include/triEle.h"
//#include "grid/include/node.h"
#include "gridGenerate/include/coord.h"
#include <unordered_map>
#include <map>


//ʹ��CGAL���ɵ������Σ������ʹ������̫��ʱ�ˡ���Ҫ�޸ġ�
// 20240604 22��29
// ����Edge�Ĺ�ϣ����
struct EdgeHash {
	std::size_t operator()(const std::shared_ptr<TriEdge>& k) const;
};

struct EdgeEqual {
	bool operator()(const std::shared_ptr<TriEdge>& lhs, const std::shared_ptr<TriEdge>& rhs) const {
		return *lhs == *rhs;
	}
};

class TriBase
{
private:
	using triElePtr_ = std::shared_ptr<TriEle>;
	using triEleMap_ = std::unordered_map<size_t, triElePtr_>;
	using cellVec_ = std::vector< std::vector<size_t> >;
	using EdgeTable = std::unordered_map<std::shared_ptr<TriEdge>, std::vector<triElePtr_>, EdgeHash, EdgeEqual>;

	std::string dividType_;
	triEleMap_ triEleVec_;
	std::vector<TriEle> notIdealTriVec_;
	cellVec_ cell_;
	EdgeTable edgeTable_;// �߱�ʹ��unordered_map���洢������Edge��ֵ�ǰ����ñߵ�����������ָ������б�
	std::unordered_map<size_t, TriNode> TriNodeMap_;
	std::vector<size_t> cgalPointsIndexVec_;
	std::pair<std::vector<size_t>, std::unordered_map<size_t, size_t>> vertex2tecplotIndex_;//first:�������; second:�ֲ�����
	//std::map<std::string, std::vector<size_t>> splitIndex_;
private:
	void removeTriFromEdgeTable(const triElePtr_& triEle_ptr);
	void removeTriEle(size_t triEleKey);
	void addTriEle(const triElePtr_& triEle);
	void addEdgeTableItem(const triElePtr_& triEle_ptr);
	void OutputConnectGrid2Tecplot(const std::vector<size_t> pointIndexVec, cellVec_& triVec,
		const std::unordered_map<size_t, size_t>& vertex2tecplotIndex, const std::string title) const;

public:
	TriBase();
	TriBase(std::vector<CoordPtr>& outterNode,  std::vector<CoordPtr>& innerNode, const std::string& typeName);
	void genCDT();
	void filterTri();
	void genEdgeTable();
	void reconstructNotIdealTriGrid();
	void reconstructNotIdealTriGrid(const TriEdge& edge, const std::vector<triElePtr_>& triPtr);
	void genTriAndQuad();
	void genTriAndQuad(const TriEdge& edge, const std::vector<triElePtr_>& triangles);
	void GenCompleteCell();   //�洢�ı��ε�Ԫ�������ε�Ԫ������
	void gencell(const triEleMap_& content);

	triElePtr_& GetTriEle(size_t index) { return triEleVec_.at(index); };
	triEleMap_& GetTriEleVec() { return triEleVec_; };
	const Coord& index2Coor(const size_t index) const { return TriNodeMap_.at(index).getcoord(); };
	std::vector<TriEle>& getNotIdealTriVec()  { return notIdealTriVec_; };
	void addNotIdealTri(const TriEle& triEle) { notIdealTriVec_.push_back(triEle); };
	const EdgeTable& getEdgeTable() const { return edgeTable_; };
	cellVec_& getCell() { return cell_; };
	//std::map<std::string, std::vector<size_t>>& getsplitIndex() { return splitIndex_; }
	void boostMesh();

};


#endif // !TRI_ELEBASE_H