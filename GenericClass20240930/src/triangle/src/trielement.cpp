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

#include "triangle/include/triEle.h"
#include <stdexcept>
#include <cmath>
#include <memory>


TriEle::TriEle(const std::vector<TriNode>& ptrVec, size_t index, TriEleTag ttag, TriEleType ttype)
	:TriIndex_(index), triEleTag_(ttag), trieleType_(ttype)
{
	size_t a = ptrVec.at(0).getindex();
	size_t b = ptrVec.at(1).getindex();
	size_t c = ptrVec.at(2).getindex();
	triVertexIndex_ = { a,b,c };
	triEdges_ = { std::make_shared<TriEdge>(a, b),std::make_shared<TriEdge>(b, c),std::make_shared<TriEdge>(c, a)};
	triNodesPtr_ = ptrVec;
	if (ttype== TriEleType::Notset)
		trieleType_ = judgeTriType();
	
}


bool TriEle::isQualifiedTri() const
{
	if (isCollinear())
		return false;
	else
		return true;
}

/*-------------------------------------------------------------------
*  Function: areCollinear(Coordinate p1, Coordinate p2, Coordinate p3)
*  Purpose: 检查三角形的三个点是否共线
*  Arguments:
*   Coordinate p1,p2,p3 - 三角形的三个顶点
*  Returns:
*    bool - 是否共线的布尔值
-------------------------------------------------------------------*/
bool TriEle::isCollinear() const
{
	//判断是否为合格的三角形
	TriNode p1 = triNodesPtr_.at(0);
	TriNode p2 = triNodesPtr_.at(1);
	TriNode p3 = triNodesPtr_.at(2);
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

//设置三角形的类型-----------------------------------
TriEleType TriEle::judgeTriType() const
{
	size_t out_num = 0;
	size_t in_num = 0;
	for (auto& perit : triNodesPtr_) {
		if (perit.gettype() == TriNodeType::outter)
			out_num++;
		else if (perit.gettype() == TriNodeType::inner)
			in_num++;
		else
		{
			//spdlog::error("curretn tri node type error!");
		}
	}
	//3个顶点都在外部边界上
	if (out_num == 3 && in_num == 0) {
		return TriEleType::Outer_only;
	}
	//2个在外部边界;1个都在内部边界
	else if (out_num == 2 && in_num == 1) {
		return TriEleType::Outer;
	}
	//1个在外部边界;2个都在内部边界
	else if (out_num == 1 && in_num == 2) {
		return TriEleType::Inner;
	}
	//3个顶点都在内部边界上
	else if (out_num == 0 && in_num == 3) {
		return TriEleType::Inner_only;
	}
	else {
		//spdlog::error("设置三角形类型时出错！！! ");
		throw std::runtime_error("Caught an exception!");
	}
}