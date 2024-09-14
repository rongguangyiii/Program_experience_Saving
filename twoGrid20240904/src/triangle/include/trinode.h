/*---------------------------------------------------------------------------
	ZaRan	-	A Totally Automatic CFD Software
	Copyright (C) ,Since 2024
-------------------------------------------------------------------------------
License
	This file is part of ZaRan.

!	@file	cdtPoint.h
!	@brief	Generate ZaRan Grid.
!	@author	Liu Guangying.
!   @date  2024.06.06
\*---------------------------------------------------------------------------*/

#ifndef TRINODE_H_
#define TRINODE_H_
//#include"Coordinate/include/Coordinate.h"
#include "gridGenerate/include/coord.h"

enum class TriNodeType
{
	Notset,
	outter,
	inner,
};

class TriNode
{
public:
	TriNode() :coord_({ 0,0,0 }), index_(0), type_(TriNodeType::Notset) {};
	TriNode(const Coord& in_coor, const size_t& id, TriNodeType ttype = TriNodeType::Notset)
		:coord_(in_coor), index_(id), type_(ttype) {};
	~TriNode() {};
	const Coord& getcoord()const { return coord_; }
	TriNodeType gettype()const { return type_; }
	const size_t& getindex() const { return index_; }
private:
	Coord coord_;
	TriNodeType type_;
	size_t index_;
};




#endif // !CDTPOINT_H_
