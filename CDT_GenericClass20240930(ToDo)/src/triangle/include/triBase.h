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
#ifndef TRI_ELEBASE_H  
#define TRI_ELEBASE_H  
#pragma once
#include "triangle/include/triEle.h"
//#include "grid/include/node.h"
#include <unordered_map>
#include <filesystem>  
#include <map>


//使用CGAL生成的三角形，这个类使用起来太耗时了。需要修改。
// 20240604 22：29
// 用于Edge的哈希函数
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
	EdgeTable edgeTable_;// 边表使用unordered_map来存储，键是Edge，值是包含该边的三角形智能指针对象列表
	std::unordered_map<size_t, TriNode> TriNodeMap_;
	std::vector<size_t> cgalPointsIndexVec_;
	std::pair<std::vector<size_t>, std::unordered_map<size_t, size_t>> vertex2tecplotIndex_;//first:点的索引; second:局部索引
	std::map<std::string, std::vector<size_t>> splitIndex_;
private:
	void removeTriFromEdgeTable(const triElePtr_& triEle_ptr);
	void removeTriEle(size_t triEleKey);
	void addTriEle(const triElePtr_& triEle);
	void addEdgeTableItem(const triElePtr_& triEle_ptr);
	void OutputConnectGrid2Tecplot(const std::vector<size_t> pointIndexVec, cellVec_& triVec,
		const std::unordered_map<size_t, size_t>& vertex2tecplotIndex, const std::string title) const;

public:
	TriBase() = delete;
	TriBase(std::vector<TriCoord>& outterNode,  std::vector<TriCoord>& innerNode, const std::string& typeName);
	void genCDT();
	void filterTri();
	void genEdgeTable();
	void reconstructNotIdealTriGrid();
	void reconstructNotIdealTriGrid(const TriEdge& edge, const std::vector<triElePtr_>& triPtr);
	void genTriAndQuad();
	void genTriAndQuad(const TriEdge& edge, const std::vector<triElePtr_>& triangles);
	void GenCompleteCell();   //存储四边形单元和三角形单元的索引
	void gencell(const triEleMap_& content);

	triElePtr_& GetTriEle(size_t index) { return triEleVec_.at(index); };
	triEleMap_& GetTriEleVec() { return triEleVec_; };
	const TriCoord& index2Coor(const size_t index) const { return TriNodeMap_.at(index).getcoord(); };
	std::vector<TriEle>& getNotIdealTriVec()  { return notIdealTriVec_; };
	void addNotIdealTri(const TriEle& triEle) { notIdealTriVec_.push_back(triEle); };
	const EdgeTable& getEdgeTable() const { return edgeTable_; };
	cellVec_& getCell() { return cell_; };
	std::map<std::string, std::vector<size_t>>& getsplitIndex() { return splitIndex_; }
	void boostMesh();

public:
	//基本方法函数
	void reorderPointsIndex(std::vector<std::pair<TriCoord, size_t>>& PointLists, const TriCoord& modelCenter);
	bool isCollinear(const TriCoord& p1, const TriCoord& p2, const TriCoord& p3);
	void createFolder(const std::string& basePath)const;
	bool isPointInPolygon_moreOuter(const TriCoord& p, const std::vector<TriCoord>& polygon);
	bool isPointInPolygon_moreInner(const TriCoord& p, const std::vector<TriCoord>& polygon);
	void reorderNodes(std::vector<TriCoord>& NodeList);
};


#endif // !TRI_ELEBASE_H