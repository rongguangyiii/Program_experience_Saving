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
*  Purpose: ��������ε��������Ƿ���
*  Arguments:
*   Coordinate p1,p2,p3 - �����ε���������
*  Returns:
*    bool - �Ƿ��ߵĲ���ֵ
-------------------------------------------------------------------*/
bool TriEle::isCollinear() const
{
	//�ж��Ƿ�Ϊ�ϸ��������
	TriNode p1 = triNodesPtr_.at(0);
	TriNode p2 = triNodesPtr_.at(1);
	TriNode p3 = triNodesPtr_.at(2);
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

//���������ε�����-----------------------------------
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
	//3�����㶼���ⲿ�߽���
	if (out_num == 3 && in_num == 0) {
		return TriEleType::Outer_only;
	}
	//2�����ⲿ�߽�;1�������ڲ��߽�
	else if (out_num == 2 && in_num == 1) {
		return TriEleType::Outer;
	}
	//1�����ⲿ�߽�;2�������ڲ��߽�
	else if (out_num == 1 && in_num == 2) {
		return TriEleType::Inner;
	}
	//3�����㶼���ڲ��߽���
	else if (out_num == 0 && in_num == 3) {
		return TriEleType::Inner_only;
	}
	else {
		//spdlog::error("��������������ʱ������! ");
		throw std::runtime_error("Caught an exception!");
	}
}