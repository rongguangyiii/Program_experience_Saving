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

enum class ConnectionType
{
	/* 1. poly: 有序凸多边形；2. line:有序直线; 3. point:无序点*/
	notset,
	point_poly,
	point_line,
	point_point,
	poly_poly,
};

class TriBase
{
	using triElePtr_ = std::shared_ptr<TriEle>;
	using triEleMap_ = std::unordered_map<size_t, triElePtr_>;
	using cellVec_ = std::vector< std::vector<size_t> >;
	using EdgeTable = std::unordered_map<std::shared_ptr<TriEdge>, std::vector<triElePtr_>, EdgeHash, EdgeEqual>;
public:
	TriBase() :connctType_(ConnectionType::notset), isoutFile(true) {};
	TriBase(const ConnectionType& type);//Debug 示例
	TriBase(std::vector<TriNode>& outterNode,  std::vector<TriNode>& innerNode, const ConnectionType& typeName);
	void genCDT();
	void filterTri();
	void genEdgeTable();
	void checkConnectiontype()const;
	void reconstructNotIdealTriGrid();
	void reconstructNotIdealTriGrid(const TriEdge& edge, const std::vector<triElePtr_>& triPtr);
	void genTriAndQuad();
	void genTriAndQuad(const TriEdge& edge, const std::vector<triElePtr_>& triangles);
	void GenCompleteCell();   //存储四边形单元和三角形单元的索引
	void gencell(const triEleMap_& content);
	void setConnectionType(const ConnectionType& type) { connctType_ = type; }
	triElePtr_& GetTriEle(size_t index) { return triEleVec_.at(index); };
	triEleMap_& GetTriEleVec() { return triEleVec_; };
	const TriNode& index2Coor(const size_t index) const { return TriNodeMap_.at(index); };
	std::vector<TriEle>& getNotIdealTriVec()  { return notIdealTriVec_; };
	void addNotIdealTri(const TriEle& triEle) { notIdealTriVec_.push_back(triEle); };
	const EdgeTable& getEdgeTable() const { return edgeTable_; };
	cellVec_& getCell() { return cell_; };
	std::map<std::string, std::vector<size_t>>& getsplitIndex() { return splitIndex_; }
	std::vector<size_t> getRemarkIdx() { return reMarkIdx_; }
	void boostMesh();
public:
	//基本方法函数
	void insertPoint(std::vector<TriNode> nodes, TriNodeType nodetype);
	void insertCtrlLine(std::vector<TriNode> nodes, TriNodeType nodetype = TriNodeType::inner);
	void insertCtrlPoly(std::vector<TriNode> nodes, TriNodeType nodetype);
	void setfileOnOff(bool onoff) { isoutFile = onoff; };
	void reorderPointsIndex(std::vector<std::pair<TriNode, size_t>>& PointLists, const TriNode& modelCenter);
	bool isCollinear(const TriNode& p1, const TriNode& p2, const TriNode& p3);
	void createFolder(const std::string& basePath)const;
	bool isPointInPolygon_moreOuter(const TriNode& p, const std::vector<TriNode>& polygon);
	bool isPointInPolygon_moreInner(const TriNode& p, const std::vector<TriNode>& polygon);
	void reorderNodes(std::vector<TriNode>& NodeList);
private:
	TriEleType judgeTriType(std::vector<TriNode> triNodes) const;
	void isRemarkBack();
	void removeTriFromEdgeTable(const triElePtr_& triEle_ptr);
	void removeTriEle(size_t triEleKey);
	void addTriEle(const triElePtr_& triEle);
	void addEdgeTableItem(const triElePtr_& triEle_ptr);
	void OutputConnectGrid2Tecplot(const std::vector<size_t> pointIndexVec, cellVec_& triVec,
		const std::unordered_map<size_t, size_t>& vertex2tecplotIndex, const std::string title) const;
private:
	bool isoutFile;
	ConnectionType connctType_;
	triEleMap_ triEleVec_;
	std::vector<TriEle> notIdealTriVec_;
	cellVec_ cell_;
	EdgeTable edgeTable_;// 边表使用unordered_map来存储，键是Edge，值是包含该边的三角形智能指针对象列表
	std::vector<size_t> reMarkIdx_;
	std::unordered_map<size_t, TriNode> TriNodeMap_;
	std::vector<size_t> cgalPointsIndexVec_;
	std::pair<std::vector<size_t>, std::unordered_map<size_t, size_t>> vertex2tecplotIndex_;//first:点的索引; second:局部索引
	std::map<std::string, std::vector<size_t>> splitIndex_;
};


#endif // !TRI_ELEBASE_H

/*----------------------------------------------------------------------------------------
* USEAGE:
* 
#include "triangle/include/triBase.h"
#include <vector>

int main() {
	//--------------------------------------------------------
	//为了构建正确的CDT，这里提供三种类型，outter-inner:
	//point_poly; point_line;  point_point;
	//完善CDT,主要增加了拓扑优化，反馈更新背景网格点的功能
	//--------------------------------------------------------
	TriBase test;
	test.setfileOnOff(true);
	test.setConnectionType(ConnectionType::point_line);
	std::vector<TriNode>poly_1_outter = { {1,2,0,1}, {3,2,0,2}, {7,2,0,3}, {7,4,0,4}, {8,7,0,5}, {5,8,0,6}, {2,6,0,7}, {3,4,0,8} };
	std::vector<TriNode>poly_2_inner = { {4,4,0,9}, {6,3,0,10}, {6,6,0,11}, {4,6,0,12}, {5,5,0,13} };
	test.insertPoint(poly_1_outter, TriNodeType::outter);
	//test.insertPoint(poly_2_inner, TriNodeType::inner);
	//test.insertCtrlPoly(poly_1_outter,TriNodeType::outter);
	test.insertCtrlLine(poly_2_inner,TriNodeType::inner);
	//test.insertCtrlPoly(poly_2_inner, TriNodeType::inner);
	test.boostMesh();

	return 0;

}
----------------------------------------------------------------------------------------*/