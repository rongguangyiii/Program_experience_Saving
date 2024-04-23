#pragma once
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <cmath>
#include "triInfo.h"

TriEdge::TriEdge(int v1, int v2)
{
	// ȷ���ߵĶ���˳���������
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
	// ʹ��FNV-1a��ϣ�㷨�ı��壬��������������ͨ����������  
	std::size_t hash = 2166136261; // FNV offset basis  
	hash ^= k.v1;
	hash *= 16777619; // FNV prime  
	hash ^= k.v2;
	hash *= 16777619; // �ٴλ��  
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

// �жϵ��Ƿ��ڶ�����ڲ��������߽�Ͷ��㣩  
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
        // �����Ƿ�պ��붥���غ�  
        if ((std::fabs(p.x - polygon[i].x) < std::numeric_limits<double>::epsilon()) &&
            (std::fabs(p.y - polygon[i].y) < std::numeric_limits<double>::epsilon())) {
            return true; // ���붥���غϣ���Ϊ���ڶ�����ڲ�  
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
	// ���������ε�������
	std::vector<TriEdge> edges{
		TriEdge(triedge_index.at(0), triedge_index.at(1)),
		TriEdge(triedge_index.at(1), triedge_index.at(2)),
		TriEdge(triedge_index.at(2), triedge_index.at(0))
	};

	// ��ÿ���ߣ�������ӵ��߱���
	for (const TriEdge& edge : edges) {
		edgeTable_[edge].push_back(triEle);
	}

}

std::vector<size_t> TriBase::returnTecplotIndex()
{
	return tecplotIndex_;
}
