/*---------------------------------------------------------------------------
	ZaRan	-	A Totally Automatic CFD Software
	Copyright (C) ,Since 2020
-------------------------------------------------------------------------------
License
	This file is part of ZaRan.

!	@file	tri_info.h
!	@brief	Generate ZaRan Grid.
!	@author	Liu Guangying.
!   @date  2024.04.20
\*---------------------------------------------------------------------------*/
#ifndef TRI_ELE_H  
#define TRI_ELE_H 
#pragma once
//#include"Coordinate/include/Coordinate.h"
#include "gridGenerate/include/coord.h"
#include"triangle/include/triedge.h"
#include"triangle/include/trinode.h"
#include<memory>



//class triBase;


enum class TriEleTag
{
	unuse,
	used,
};

enum class TriEleType
{
	Notset,
	Outer_only,
	Outer,
	Inner,
	Inner_only,

};


class TriEle 
{
private:
	size_t TriIndex_;
	TriEleType trieleType_;
	TriEleTag triEleTag_;
	std::vector<std::shared_ptr<TriEdge>> triEdges_;
	Coord triCenter_;
	std::vector<size_t> triVertexIndex_;
	std::vector<TriNode> triNodesPtr_;

public:
	TriEle() = delete;
	TriEle(const std::vector<TriNode>&, size_t index, 
		TriEleTag ttag = TriEleTag::unuse, TriEleType ttype = TriEleType::Notset);
	// 拷贝构造函数
	TriEle(const TriEle& other)
		: TriIndex_(other.TriIndex_), trieleType_(other.trieleType_),
		triEleTag_(other.triEleTag_), triEdges_(other.triEdges_),
		triCenter_(other.triCenter_), triVertexIndex_(other.triVertexIndex_),
		triNodesPtr_(other.triNodesPtr_) {}

	// 移动构造函数
	TriEle(TriEle&& other) noexcept
		: TriIndex_(std::move(other.TriIndex_)), trieleType_(std::move(other.trieleType_)),
		triEleTag_(std::move(other.triEleTag_)), triEdges_(std::move(other.triEdges_)),
		triCenter_(std::move(other.triCenter_)), triVertexIndex_(std::move(other.triVertexIndex_)),
		triNodesPtr_(std::move(other.triNodesPtr_)) {}

	const std::vector<size_t>& GetTriEleVertex() const { return triVertexIndex_; };
	size_t GetTriEleIndex() const { return TriIndex_; };
	void SetTriEleIndex(size_t currentTriIndex) { TriIndex_ = currentTriIndex; };
	const Coord& getTriCenter() const { return triCenter_; };
	void setTriCenter(const Coord& baryCenter) { triCenter_ = baryCenter; };
	TriEleTag getTriEleTag() const { return triEleTag_; };
	TriEleType getTriEleType() const { return trieleType_; };
	void setTriEleTag(TriEleTag triEleTag) { triEleTag_ = triEleTag; };
	void setTriEleType(TriEleType trieleType) { trieleType_ = trieleType; };
	std::vector<std::shared_ptr<TriEdge>>& getTriEdgesRef() { return triEdges_; };
	TriEleType judgeTriType()const;
	bool isQualifiedTri() const;
	bool isCollinear() const;

};


#endif // TRI_ELE_H



