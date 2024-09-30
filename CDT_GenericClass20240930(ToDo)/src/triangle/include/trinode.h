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
#ifndef TRINODE_H_
#define TRINODE_H_
#include"triangle/include/triCoord.h"


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
	TriNode(const TriCoord& in_coor, const size_t id, TriNodeType ttype = TriNodeType::Notset)
		:coord_(in_coor), index_(id), type_(ttype) {};
	~TriNode() {};
	const TriCoord& getcoord()const { return coord_; }
	TriNodeType gettype()const { return type_; }
	const size_t getindex() const { return index_; }
private:
	TriCoord coord_;
	TriNodeType type_;
	size_t index_;
};




#endif // !CDTPOINT_H_
