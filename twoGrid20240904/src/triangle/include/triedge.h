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

// 定义边结构，使用顶点索引表示
class TriEdge {
private:

public:
	size_t vA, vB;
	//实际上这里还应该有TriNode信息,但并没有提供额外有用的信息，就没有写。
public:
	TriEdge() :vA(0), vB(0) {};
	TriEdge(size_t vp1, size_t vp2) :vA(std::min(vp1, vp2)), vB(std::max(vp1, vp2)) {}
	// 为了使用Edge作为unordered_map的键，需要重载==和hash函数
	bool operator==(const TriEdge& other) const {
		return vA == other.vA && vB == other.vB;
	}

};



#endif // !TRI_EDGE_H
