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
#ifndef TRI_EDGE_H
#define TRI_EDGE_H
//#include<cmath>
#include<algorithm>
#include<memory>

// ����߽ṹ��ʹ�ö���������ʾ
class TriEdge {
private:

public:
	size_t vA, vB;
	//ʵ�������ﻹӦ����TriNode��Ϣ,����û���ṩ�������õ���Ϣ����û��д��
public:
	TriEdge() :vA(0), vB(0) {};
	TriEdge(size_t vp1, size_t vp2) :vA(std::min(vp1, vp2)), vB(std::max(vp1, vp2)) {}
	// Ϊ��ʹ��Edge��Ϊunordered_map�ļ�����Ҫ����==��hash����
	bool operator==(const TriEdge& other) const {
		return vA == other.vA && vB == other.vB;
	}

};



#endif // !TRI_EDGE_H
