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
	edgeType_ = EdgeType::Notset;
}
TriEdge::TriEdge()
{
	v1 = 0;
	v2 = 0;
	edgeType_ = EdgeType::Notset;
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

void TriBase::updateTriEleIndex()
{
	size_t currentTriIndex = 0;
	for (auto& triEle : triEleVec_)
	{
		triEle.reSetTriEleIndex(currentTriIndex++);
	}

	//只要更新过triEle中的相关信息，就需要更新边表
	updateEdgeTable();

	return;
}


TriEle::TriEle()
{
	a_ = 0;
	b_ = 0;
	c_ = 0;
	EleIndex_ = 0;
	trieleType_ = TriEleType::Notset;
	triEleTag_ = TriEleTag::unuse;
}
TriEle::TriEle(size_t a, size_t b, size_t c, size_t index)
{
	a_ = a;
	b_ = b;
	c_ = c;
	EleIndex_ = index;
	trieleType_ = TriEleType::Notset;
	triEleTag_ = TriEleTag::unuse;
	std::vector<TriEdge> triEdges{
		TriEdge(a_, b_),
		TriEdge(b_, c_),
		TriEdge(c_, a_)
	};
	triEdges_ = triEdges;
}
std::vector<size_t> TriEle::GetTriEleVertex() const
{
	std::vector<size_t> TriVertexVec;
	TriVertexVec.push_back(a_);
	TriVertexVec.push_back(b_);
	TriVertexVec.push_back(c_);
	return TriVertexVec;
}
void TriEle::setTriEleTag(TriEleTag triEleTag)
{
	triEleTag_ = triEleTag;
}

TriEleTag TriEle::getTriEleTag() const
{
	return triEleTag_;
}

TriEleType TriEle::getTriEleType() const
{
	return trieleType_;
}

std::vector<TriEdge>& TriEle::getTriEdges()
{
	return triEdges_;
}

void TriEle::setTriEleType(TriEleType trieleType)
{
	trieleType_ = trieleType;
}
// 判断点是否在多边形内部（包括边界和顶点）  
bool isPointInPolygon(const Point& p, const std::vector<Point>& polygon) {
	bool inside = false;
	size_t j = polygon.size() - 1;
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
	// 对每个边，将其添加到边表中
	for (const TriEdge& edge : triEle.getTriEdges()) {
		edgeTable_[edge].push_back(triEle);
	}

}

void TriBase::updateEdgeTable()
{
	edgeTable_.clear();
	for (auto& triEle : triEleVec_)
	{
		addTri2EdgeTable(triEle);
	}
	return;
}

std::vector<size_t> TriBase::returnTecplotIndex()
{
	return tecplotIndex_;
}

const TriBase::EdgeTable& TriBase::getEdgeTable() const
{
	return edgeTable_;
}

void TriBase::GenTriAndquad()
{
	auto& triVec = GetTriEleVec();
	for (size_t iTri = 0; iTri < triVec.size(); ++iTri)
	{
		auto& cur_tri = triVec.at(iTri);
		TriEleTag cur_triTag = cur_tri.getTriEleTag();
		//如果当前三角形已经用于连结过，就直接进行下一次循环。
		if (TriEleTag::used == cur_triTag)
			continue;

		bool isReashape = false;
		size_t cur_triIndex = cur_tri.GetTriEleIndex();
		TriEleType cur_triType = cur_tri.getTriEleType();
		for (TriEdge& edge : cur_tri.getTriEdges())
		{
			if (EdgeType::unique == edge.getEdgeType())
				continue;
			if (isReashape)
				continue;

			/*auto& edgetable = getEdgeTable();
			auto& it = edgetable.find(edge);*/
			auto it = getEdgeTable().find(edge);
			const std::vector<TriEle>& triangles = it->second;
			//对于当前三角形的内部边，检查相邻的两个三角形是否可以重构为四边形
			for (auto& var : triangles)
			{
				if (cur_triIndex == var.GetTriEleIndex())
					continue;
				if (cur_triType == var.getTriEleType())
					continue;
				if (TriEleTag::used == var.getTriEleTag())
					continue;
				/*auto target_Index = var.GetTriEleIndex();
				auto target_type = var.getTriEleType();
				auto target_tag = triVec.at(var.GetTriEleIndex()).getTriEleTag();*/
				isReashape = true;
			}

			if (isReashape)
			{
				//如果当前三角形有符合要求的相邻三角形，就进行重构为四边形，并且同时标记相邻三角形为TriEleTag::used。
				GenTriAndquad(edge, triangles);
				//只要更新过triEle的相关信息，就需要更新边表,这里主要更新相邻三角形为TriEleTag::used
				updateEdgeTable();
				break;
			}
		}//end for edge
	}//end for iTri

	//所有三角形循环结束后，再处理仍然未被使用的三角形
	GenCompleteCell();

	return;
}

//判断三角形的三个点其中两个在外环上还是内环上，并设置三角形的类型
void TriBase::updateTriGridType()
{
	//step 1:获取背景网格和非结构网格的洞边界节点(在非结构网格下)索引
	std::vector<size_t> poly_outer_index = { 0, 1, 2, 3, 4, 5, 6, 7, };
	size_t num_1 = 0;
	size_t num_2 = 0;
	std::vector<TriEle>& TriVec = GetTriEleVec();
	//step 2: 检查三角形的三个顶点的其中两个顶点在外环上还是内环上
	for (auto& var : TriVec)
	{
		auto triVertex = var.GetTriEleVertex();

		size_t have_count = 0;
		size_t notHave_count = 0;
		//每个三角形的三个顶点循环
		for (auto perVertex : triVertex)
		{
			// 检查是否找到了元素  
			auto it = std::find(poly_outer_index.begin(), poly_outer_index.end(), perVertex);
			if (it != poly_outer_index.end())
				have_count++;
			else
				notHave_count++;
		}
		//TriEle的TrieleType_默认为Notset.
		TriEleType currentTriType = TriEleType::Notset;
		if (have_count == 2 && notHave_count == 1 && have_count + notHave_count == triVertex.size())//理想三角形
		{
			//当前三角形的三个顶点中有两个在外环上，一个在内环上
			currentTriType = TriEleType::Outer;
			var.setTriEleType(currentTriType);
			num_1++;
		}
		else if (have_count == 1 && notHave_count == 2 && have_count + notHave_count == triVertex.size())//三个顶点都在外部边界上或都在内部边界上
		{
			//当前三角形的三个顶点中有一个在外环上，两个在内环上
			currentTriType = TriEleType::Inner;
			var.setTriEleType(currentTriType);
			num_2++;
		}
		else
		{
			std::cout << "设置三角类型出错！！！ " << std::endl;
			std::cout << "程序已退出，退出错误编号：335" << std::endl;
			exit(335);
		}
	}

	//只要更新过triEle中的相关信息，就需要更新边表
	updateEdgeTable();

	return;
}

void TriBase::GenTriAndquad(const TriEdge& edge, const std::vector<TriEle>& triangles)
{
	std::vector<size_t> tri_1_vertex = triangles.at(0).GetTriEleVertex();
	std::vector<size_t> tri_2_vertex = triangles.at(1).GetTriEleVertex();
	size_t coindex_a = edge.v1;
	size_t coindex_b = edge.v2;
	size_t tri_1_diff = 0;
	size_t tri_2_diff = 0;
	//找到第三个和coindex_a，coindex_b不同的顶点
	for (auto it : tri_1_vertex)
	{
		if (it != coindex_a && it != coindex_b)
			tri_1_diff = it;
	}
	for (auto it : tri_2_vertex)
	{
		if (it != coindex_a && it != coindex_b)
			tri_2_diff = it;
	}

	std::vector<size_t> quad{ tri_1_diff, tri_2_diff, coindex_a, coindex_b };
	std::vector<std::pair<Point, size_t >> quad_pair;
	Point average;
	for (auto it : quad)
	{
		Point temcoor = all_points_.at(it);
		quad_pair.push_back(std::make_pair(all_points_.at(it), it));
		average += (temcoor / quad.size());
	}

	reorderPointsIndex(quad_pair, average);//按逆时针重新排序
	quad.clear();
	for (auto& it : quad_pair)
		quad.push_back(it.second);

	//---------------向非结构网格中添加新单元------------------------------------------------------------
	cell_.push_back(quad);
	//------------------------------------------------------------------------------------------------

	//存在的问题：1.边表对应的三角形没有设置过类型TriEleType，可能是未用引用，
	//          2.已经使用的三角形，未进行TriEleTag::used标记。怎么标记？
	auto& triVec = GetTriEleVec();
	for (auto& var : triangles)
	{
		if (triVec.at(var.GetTriEleIndex()).getTriEleTag() == TriEleTag::used)
		{
			std::cout << "索引为( " << var.GetTriEleIndex() << " )的三角形转换四边形出错！！！" << std::endl;
			std::cout << "程序已退出，退出错误编号：336" << std::endl;
			exit(336);
		}
		triVec.at(var.GetTriEleIndex()).setTriEleTag(TriEleTag::used);
	}

	return;
}

void TriBase::GenCompleteCell()
{
	for (auto& triEle : triEleVec_)
	{
		if (TriEleTag::used == triEle.getTriEleTag())
			continue;
		cell_.push_back(triEle.GetTriEleVertex());
	}

	return;
}

void TriBase::calTriBaryCenter()
{
	for (auto& triEle : triEleVec_)
	{
		auto triVertex = triEle.GetTriEleVertex();
		//计算三角形的重心
		Point ave_triCenter;
		double x_ave = 0;
		double y_ave = 0;
		for (auto perVertex : triVertex)
		{
			x_ave += all_points_.at(perVertex).x / triVertex.size();
			y_ave += all_points_.at(perVertex).y / triVertex.size();
		}
		ave_triCenter = { x_ave, y_ave };
		triEle.setBaryCenter(ave_triCenter);
	}
	return;
}

void TriBase::setBaryCenter(TriEle& tri)
{
	auto triVertex = tri.GetTriEleVertex();
	//计算三角形的重心
	Point ave_triCenter;
	double x_ave = 0;
	double y_ave = 0;
	for (auto perVertex : triVertex)
	{
		x_ave += all_points_.at(perVertex).x / triVertex.size();
		y_ave += all_points_.at(perVertex).y / triVertex.size();
	}
	ave_triCenter = { x_ave, y_ave };
	tri.setBaryCenter(ave_triCenter);
}



void TriBase::WriteTecplotFile(const std::string& title)
{
	std::ofstream outfile(title + "_grid" + ".dat");

	// 写入文件头
	outfile << "TITLE = \"TriBase Tecplot Output\"" << std::endl;
	outfile << "VARIABLES = \"X\", \"Y\"" << std::endl;

	// 写入节点信息
	outfile << "ZONE N=" << all_points_.size() << ", E=" << cell_.size() << ", F=FEPOINT, ET=QUADRILATERAL" << std::endl;
	for (const auto& point : all_points_) {
		outfile << point.x << " " << point.y << std::endl;
	}

	// 写入单元信息
	for (const auto& cell : cell_) {
		if (cell.size() == 3) { // 三角形单元
			outfile << cell[0] + 1 << " " << cell[1] + 1 << " " << cell[2] + 1 << " " << cell[0] + 1 << std::endl;
		}
		else if (cell.size() == 4) { // 四边形单元
			outfile << cell[0] + 1 << " " << cell[1] + 1 << " " << cell[2] + 1 << " " << cell[3] + 1 << std::endl;
		}
	}

	outfile.close();

	return;
}
