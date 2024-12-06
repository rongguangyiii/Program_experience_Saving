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

enum class DataType
{
	/*ring: ����͹����Σ�line:����ֱ��; point:�����*/
	NUset,
	ring_ring,
	point_line,
	point_ring,
	point_point,
};

class TriBase
{
	using triElePtr_ = std::shared_ptr<TriEle>;
	using triEleMap_ = std::unordered_map<size_t, triElePtr_>;
	using cellVec_ = std::vector< std::vector<size_t> >;
	using EdgeTable = std::unordered_map<std::shared_ptr<TriEdge>, std::vector<triElePtr_>, EdgeHash, EdgeEqual>;
public:
	TriBase() :dividType_(DataType::NUset), isoutFile(true) {};
	TriBase(const DataType& type);//Debug ʾ��
	TriBase(std::vector<TriNode>& outterNode,  std::vector<TriNode>& innerNode, const DataType& typeName);
	void genCDT();
	void filterTri();
	void genEdgeTable();
	void reconstructNotIdealTriGrid();
	void reconstructNotIdealTriGrid(const TriEdge& edge, const std::vector<triElePtr_>& triPtr);
	void genTriAndQuad();
	void genTriAndQuad(const TriEdge& edge, const std::vector<triElePtr_>& triangles);
	void GenCompleteCell();   //�洢�ı��ε�Ԫ�������ε�Ԫ������
	void gencell(const triEleMap_& content);
	void setConstrainType(const DataType& type) { dividType_ = type; }
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
	//������������
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
	DataType dividType_;
	triEleMap_ triEleVec_;
	std::vector<TriEle> notIdealTriVec_;
	cellVec_ cell_;
	EdgeTable edgeTable_;// �߱�ʹ��unordered_map���洢������Edge��ֵ�ǰ����ñߵ�����������ָ������б�
	std::vector<size_t> reMarkIdx_;
	std::unordered_map<size_t, TriNode> TriNodeMap_;
	std::vector<size_t> cgalPointsIndexVec_;
	std::pair<std::vector<size_t>, std::unordered_map<size_t, size_t>> vertex2tecplotIndex_;//first:�������; second:�ֲ�����
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
	//Ϊ�˹�����ȷ��CDT�������ṩ�������ͣ�outter-inner:
	//ring_ring;  point_point; point_line; point_ring;
	//����CDT,��Ҫ�����������Ż����������±��������Ĺ���
	//--------------------------------------------------------
TriBase test;
test.setConstrainType(DataType::ring_ring);
std::vector<TriNode>poly_1_outter = { {1,2,0,1}, {3,2,0,2}, {7,2,0,3}, {7,4,0,4}, {8,7,0,5}, {5,8,0,6}, {2,6,0,7}, {3,4,0,8} };
std::vector<TriNode>poly_2_inner = { {4,4,0,9}, {6,3,0,10}, {6,6,0,11}, {4,6,0,12}, {5,5,0,13} };
//test.insertPoint(poly_1_outter, TriNodeType::outter);
//test.insertPoint(poly_2_inner, TriNodeType::inner);
test.insertCtrlPoly(poly_1_outter, TriNodeType::outter);
test.insertCtrlPoly(poly_2_inner, TriNodeType::inner);
test.setfileOnOff(true);
test.boostMesh();

return 0;

}
----------------------------------------------------------------------------------------*/