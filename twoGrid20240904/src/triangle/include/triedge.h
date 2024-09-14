/*---------------------------------------------------------------------------
	ZaRan	-	A Totally Automatic CFD Software
	Copyright (C) ,Since 2024
-------------------------------------------------------------------------------
License
	This file is part of ZaRan.

!	@file	triedge.h
!	@brief	Generate ZaRan Grid.
!	@author	Liu Guangying.
!   @date  2024.06.06
\*---------------------------------------------------------------------------*/

#ifndef TRI_EDGE_H
#define TRI_EDGE_H
#include <algorithm>
#include <vector>
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
