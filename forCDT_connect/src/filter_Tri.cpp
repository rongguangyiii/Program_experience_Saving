#pragma once
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <cmath>
#include "triInfo.h"

TriEdge::TriEdge(int v1, int v2)
{
	// 确保边的顶点顺序是排序的
	if (v1 > v2)
	{
		this->v1 = v2;
		this->v2 = v1;
	}
	else
	{
		this->v1 = v1;
		this->v2 = v2;
	}
	edgeType_ = EdgeType::unset;
}
TriEdge::TriEdge()
{
	v1 = 0;
	v2 = 0;
	edgeType_ = EdgeType::unset;
};

std::size_t EdgeHash::operator()(const TriEdge& k) const {
	// 使用FNV-1a哈希算法的变体，它对于整数输入通常表现良好  
	std::size_t hash = 2166136261; // FNV offset basis  
	hash ^= k.v1;
	hash *= 16777619; // FNV prime  
	hash ^= k.v2;
	hash *= 16777619; // 再次混合  
	return hash;
}

TriEle TriBase::GetTriEle(size_t index) const
{
	return triEleVec_.at(index);
}

void TriBase::AddTriEle(TriEle triEle)
{
	triEleVec_.push_back(triEle);
}

void TriBase::reSetAllTriEleIndex( )
{
	size_t currentTriIndex= 0;
	for (auto& triEle : triEleVec_)
	{
		triEle.reSetTriEleIndex(currentTriIndex++);
	}
}


TriEle::TriEle()
{
	a_ = 0;
	b_ = 0;
	c_ = 0;
	EleIndex_ = 0;
}
TriEle::TriEle(size_t a, size_t b, size_t c, size_t index)
{
	a_ = a;
	b_ = b;
	c_ = c;
	EleIndex_ = index;
}
std::vector<size_t> TriEle::GetTriEleVertex() const
{
	std::vector<size_t> TriVertexVec;
	TriVertexVec.push_back(a_);
	TriVertexVec.push_back(b_);
	TriVertexVec.push_back(c_);
	return TriVertexVec;
}

// 判断点是否在多边形内部（包括边界和顶点）  
bool isPointInPolygon(const Point& p, const std::vector<Point>& polygon) {
    bool inside = false;
    int j = polygon.size() - 1;
    for (size_t i = 0; i < polygon.size(); i++) {
        if (((polygon[i].y <= p.y) && (p.y < polygon[j].y)) ||
            ((polygon[j].y <= p.y) && (p.y < polygon[i].y))) {
            if (polygon[i].x + (p.y - polygon[i].y) * (polygon[j].x - polygon[i].x) / (polygon[j].y - polygon[i].y) < p.x) {
                inside = !inside;
            }
        }
        // 检查点是否刚好与顶点重合  
        if ((std::fabs(p.x - polygon[i].x) < std::numeric_limits<double>::epsilon()) &&
            (std::fabs(p.y - polygon[i].y) < std::numeric_limits<double>::epsilon())) {
            return true; // 点与顶点重合，认为点在多边形内部  
        }
        j = i;
    }
    return inside;
}
void TriBase::addNotIdealTri(TriEle triEle)
{
	notIdealTriVec_.push_back(triEle);
}

void TriBase::addTri2EdgeTable(TriEle& triEle)
{
	auto&triedge_index = triEle.GetTriEleVertex();
	// 创建三角形的三个边
	std::vector<TriEdge> edges{
		TriEdge(triedge_index.at(0), triedge_index.at(1)),
		TriEdge(triedge_index.at(1), triedge_index.at(2)),
		TriEdge(triedge_index.at(2), triedge_index.at(0))
	};

	// 对每个边，将其添加到边表中
	for (const TriEdge& edge : edges) {
		edgeTable_[edge].push_back(triEle);
	}

}

std::vector<size_t> TriBase::returnTecplotIndex()
{
	return tecplotIndex_;
}
