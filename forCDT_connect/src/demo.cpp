#include <iostream>
#include <fstream> // �����ļ�����/���
#include <map>
#include <algorithm>
#include <vector>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include "triInfo.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef CGAL::Triangulation_vertex_base_with_info_2<int, K> Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;
typedef CDT::Vertex_handle Vertex_handle;
typedef CDT::Finite_faces_iterator Finite_faces_iterator;

int main() {
	TriBase triBase;
	triBase.genCDT();
	triBase.filterTriGrid();
	triBase.checkTriGrid();
	triBase.reconstructNotIdealTriGrid();
	triBase.reSetAllTriEleIndex();

	triBase.OutputTriGrid2Tecplot(triBase.returnTecplotIndex(), triBase.GetTriEleVec(), "Tri_reconstruct");
	return 0;
}


int TriBase::genCDT() {
	std::vector<Point_2> poly_1_cdt = { {1,2}, {3,2}, {7,2}, {7,4}, {8,7}, {5,8}, {2,6}, {3,4} };
	std::vector<Point_2> poly_2_cdt = { {4,4}, {6,3}, {6,6}, {4,6}, {5,5} };

	// �ϲ������㼯  
	std::vector<Point_2>all_points = poly_1_cdt;
	all_points.insert(all_points.end(), poly_2_cdt.begin(), poly_2_cdt.end());
	all_points_ = { {1,2}, {3,2}, {7,2}, {7,4}, {8,7}, {5,8}, {2,6}, {3,4},{4,4}, {6,3}, {6,6}, {4,6}, {5,5} };
	poly_outer = { {1,2}, {3,2}, {7,2}, {7,4}, {8,7}, {5,8}, {2,6}, {3,4} };
	poly_inner = { {4,4}, {6,3}, {6,6}, {4,6}, {5,5} };
	// ����һ��Լ�������ʷ�
	CDT cdt;
	std::map<Vertex_handle, int> vertex_to_index;
	int index = 0;
	for (const Point_2& p : all_points)
	{
		tecplotIndex_.push_back(index);
		Vertex_handle vh = cdt.insert(p);
		vertex_to_index[vh] = index++;
	}

	for (size_t i = 0; i < poly_1_cdt.size(); ++i)
	{
		cdt.insert_constraint(poly_1_cdt[i], poly_1_cdt[(i + 1) % poly_1_cdt.size()]);
	}
	for (size_t i = 0; i < poly_2_cdt.size(); ++i)
	{
		cdt.insert_constraint(poly_2_cdt[i], poly_2_cdt[(i + 1) % poly_2_cdt.size()]);
	}

	// ����һ�� Tecplot �ļ�
	std::ofstream tecplotFile("Tri_ini_grid.dat");
	tecplotFile << "TITLE = \"Triangle Output\"\n";
	tecplotFile << "VARIABLES = \"X\", \"Y\"\n";

	// д�붥��
	tecplotFile << "ZONE T=\"Triangles\", N=" << cdt.number_of_vertices() << ", E=" << cdt.number_of_faces() << ", F=FEPOINT, ET=TRIANGLE\n";
	for (auto vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit) {
		Point_2 p = vit->point();
		tecplotFile << p.x() << " " << p.y() << "\n";
	}

	// д��������
	size_t triIndex = 0;
	for (Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
		tecplotFile << vertex_to_index[fit->vertex(0)] + 1 << " "
			<< vertex_to_index[fit->vertex(1)] + 1 << " "
			<< vertex_to_index[fit->vertex(2)] + 1 << "\n";
		TriEle triEle(vertex_to_index[fit->vertex(0)], vertex_to_index[fit->vertex(1)], vertex_to_index[fit->vertex(2)], triIndex);
		AddTriEle(triEle);
		++triIndex;
	}

	tecplotFile.close();

	return 0;
}

void TriBase::filterTriGrid()
{

	//step 1: �����ⲿ���ڲ��ı߽�����

	std::vector<TriEle>& pickTriVec = GetTriEleVec();//���ù�ϵ��ֱ���޸�tribase_�е�triEleVec_����

	//step 2 : �ж�һ���������Ƿ����ڲ�߽������ڲ���������ڲ���ɾ��
	size_t perTri_index = 0;
	for (size_t iTri = 0; iTri < pickTriVec.size(); ++iTri)
	{
		auto triVertex = pickTriVec.at(iTri).GetTriEleVertex();
		//���������ε�����
		Point ave_triCenter;
		double x_ave = 0;
		double y_ave = 0;
		for (auto perVertex : triVertex)
		{
			x_ave += all_points_.at(perVertex).x / triVertex.size();
			y_ave += all_points_.at(perVertex).y / triVertex.size();
		}
		ave_triCenter = { x_ave, y_ave };
		pickTriVec.at(iTri).setBaryCenter(ave_triCenter);

		if (isPointInPolygon(ave_triCenter, poly_inner))
		{
			pickTriVec.erase(pickTriVec.begin() + iTri);
			--iTri; //ɾ��һ�������κ󣬺����������������ǰ�ƣ�������Ҫ��1,�����©��һ�������ε��ж�
		}
	}

	//step 3 : �ж�һ���������Ƿ������߽������ⲿ��������ⲿ��ɾ��
	for (size_t iTri = 0; iTri < pickTriVec.size(); ++iTri)
	{
		auto triVertex = pickTriVec.at(iTri).GetTriEleVertex();
		//���������ε�����
		Point ave_triCenter;
		double x_ave = 0;
		double y_ave = 0;
		for (auto perVertex : triVertex)
		{
			x_ave += all_points_.at(perVertex).x / triVertex.size();
			y_ave += all_points_.at(perVertex).y / triVertex.size();
		}
		ave_triCenter = { x_ave, y_ave };

		if (!isPointInPolygon(ave_triCenter, poly_outer))
		{
			pickTriVec.erase(pickTriVec.begin() + iTri);
			--iTri;
		}
	}

	OutputTriGrid2Tecplot(tecplotIndex_, GetTriEleVec(), "Tri_filterd");

	return;
}

void TriBase::OutputTriGrid2Tecplot(const std::vector<size_t> pointIndexVec, std::vector<TriEle> triVec, const std::string title)
{
	// ����һ�� Tecplot �ļ�
	std::ofstream tecplotFile(title + "_grid" + ".dat");
	//std::ofstream tecplotFile("output.dat");
	tecplotFile << "TITLE = \"Triangle Output\"\n";
	tecplotFile << "VARIABLES = \"X\", \"Y\"\n";

	// д�붥������
	tecplotFile << "ZONE T=\"Triangles\", N=" << pointIndexVec.size()
		<< ", E=" << triVec.size() << ", F=FEPOINT, ET=TRIANGLE\n";
	for (const auto& vit : pointIndexVec)
	{
		const auto& p = all_points_.at(vit);
		tecplotFile << p.x << " " << p.y << "\n";
	}

	// д�������ζ�������
	for (const auto& fit : triVec)
	{
		const auto& perit = fit.GetTriEleVertex();
		tecplotFile << perit.at(0) + 1 << " "
			<< perit.at(1) + 1 << " "
			<< perit.at(2) + 1 << "\n";
	}
	tecplotFile.close();

	return;
}

void TriBase::checkTriGrid()
{
	reSetAllTriEleIndex();//�������������ε�Ԫ������
	std::vector<size_t> poly_outer_bound{ 0,1,2,3,4,5,6,7 };
	std::vector<TriEle>  notIdealTriVec;
	size_t Tri_count = 0;
	for (size_t iTri = 0; iTri < triEleVec_.size(); ++iTri)
	{
		size_t have_count = 0;
		size_t notHave_count = 0;
		auto triVertex = triEleVec_.at(iTri).GetTriEleVertex();
		addTri2EdgeTable(triEleVec_.at(iTri));
		//ÿ�������ε���������ѭ��
		for (auto perVertex : triVertex)
		{
			// ����Ƿ��ҵ���Ԫ��  
			auto it = std::find(poly_outer_bound.begin(), poly_outer_bound.end(), perVertex);
			if (it != poly_outer_bound.end())
				have_count++;
			else
				notHave_count++;
		}

		if (have_count != 0 && notHave_count != 0 && have_count + notHave_count == triVertex.size())//����������
		{
			++Tri_count;
			continue;
		}
		else if (have_count == triVertex.size())//�������㶼���ⲿ�߽���
		{
			notIdealTriVec.push_back(triEleVec_.at(iTri));
			addNotIdealTri(triEleVec_.at(iTri));
		}
		else if (notHave_count == triVertex.size())//�������㶼û�����ⲿ�߽���
		{
			notIdealTriVec.push_back(triEleVec_.at(iTri));
			addNotIdealTri(triEleVec_.at(iTri));
		}
		else
		{
			std::cout << "���������ʱ�������� " << std::endl;
		}
	}

	//��һ����������Ҫ���������

	return;

}

void TriBase::reconstructNotIdealTriGrid()
{
	std::vector<TriEle> notIdealTriVec = getNotIdealTriVec();
	if (0 == notIdealTriVec.size())
		return;

	reSetAllTriEleIndex();//�������������ε�Ԫ������
	size_t Tri_count = 0;
	std::vector<size_t> toDeleteTriIndex;
	for (size_t inot = 0; inot < notIdealTriVec.size(); inot++)
	{
		auto& currentTri = notIdealTriVec.at(inot);
		auto triVertex = currentTri.GetTriEleVertex();
		std::vector<TriEdge> edges{
			TriEdge(triVertex.at(0), triVertex.at(1)),
			TriEdge(triVertex.at(1), triVertex.at(2)),
			TriEdge(triVertex.at(2), triVertex.at(0))
		};
		size_t numTri = 0;
		for (const TriEdge& edge : edges)
		{
			auto it = edgeTable_.find(edge);
			if (it != edgeTable_.end())
				numTri = it->second.size();
			if (2 == numTri)
			{
				const std::vector<TriEle>& triangles = it->second;
				//Ϊ��ֹ��ɾ�����Σ��Ƚ�Ҫɾ���������������洢����
				for each (auto var in triangles)
					toDeleteTriIndex.push_back(var.GetTriEleIndex());
				reconstructNotIdealTriGrid(edge, triangles);

			}
		}
	}

	//ɾ������������
	DeleteTri(toDeleteTriIndex);

	return;
}


void reorderPointsIndex(std::vector<std::pair<Point, size_t>>& PointLists, const Point& modelCenter)
{
	std::size_t wallNodeNum = PointLists.size();
	std::vector<std::pair<double, size_t>> angleAndIndex;
	angleAndIndex.resize(wallNodeNum);
	for (size_t iNode = 0; iNode < wallNodeNum; ++iNode)
	{
		auto& perNode = PointLists.at(iNode).first;
		Point vec = perNode - modelCenter;
		double angle = atan2(vec.y, vec.x);
		angleAndIndex[iNode] = std::make_pair(angle, iNode);
	}
	std::sort(angleAndIndex.begin(), angleAndIndex.end(),
		[](const std::pair<double, size_t>& a, const std::pair<double, size_t>& b) {return a.first < b.first; });

	std::vector<size_t> reorderedIndices;
	std::vector<Point> NodeVec_temp;
	std::vector<std::pair<Point, size_t>> NodePairVec_temp2;
	NodeVec_temp.resize(wallNodeNum);
	reorderedIndices.resize(wallNodeNum);
	for (size_t i = 0; i < wallNodeNum; ++i)
	{
		NodeVec_temp[i] = PointLists[angleAndIndex[i].second].first;
		reorderedIndices[i] = PointLists[angleAndIndex[i].second].second;
		NodePairVec_temp2.push_back(std::make_pair(NodeVec_temp[i], reorderedIndices[i]));
	}
	PointLists = NodePairVec_temp2;//���޸�ԭ�����˳��;
}

void TriBase::reconstructNotIdealTriGrid(const TriEdge& edge, const std::vector<TriEle>& triangles)
{
	std::vector<size_t> tri_1_vertex = triangles.at(0).GetTriEleVertex();
	std::vector<size_t> tri_2_vertex = triangles.at(1).GetTriEleVertex();
	size_t coindex_a = edge.v1;
	size_t coindex_b = edge.v2;
	size_t tri_1_diff = 0;
	size_t tri_2_diff = 0;
	//�ҵ���������coindex_a��coindex_b��ͬ�Ķ���
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
	std::vector<size_t> newTri_1{ tri_1_diff, tri_2_diff, coindex_a };
	std::vector<size_t> newTri_2{ tri_1_diff, tri_2_diff, coindex_b };
	std::vector<std::pair<Point, size_t >> newTri_1_pair;
	std::vector<std::pair<Point, size_t >> newTri_2_pair;
	for (auto it : newTri_1)
	{
		Point temcoor = all_points_.at(it);
		newTri_1_pair.push_back(std::make_pair(temcoor, it));
	}
	for (auto it : newTri_2)
	{
		Point temcoor = all_points_.at(it);
		newTri_2_pair.push_back(std::make_pair(temcoor, it));
	}
	reorderPointsIndex(newTri_1_pair, triangles.at(0).getBaryCenter());//����ʱ����������
	reorderPointsIndex(newTri_2_pair, triangles.at(1).getBaryCenter());
	newTri_1.clear();
	newTri_2.clear();
	for (auto& it : newTri_1_pair)
		newTri_1.push_back(it.second);
	for (auto& it : newTri_2_pair)
		newTri_2.push_back(it.second);
	TriEle triEle_1(newTri_1.at(0), newTri_1.at(1), newTri_1.at(2), GetTriEleVec().size());
	AddTriEle(triEle_1);
	TriEle triEle_2(newTri_2.at(0), newTri_2.at(1), newTri_2.at(2), GetTriEleVec().size());
	AddTriEle(triEle_2);

	return;
}

void TriBase::DeleteTri(std::vector<size_t> toDeleteTriIndex)
{
	std::sort(toDeleteTriIndex.begin(), toDeleteTriIndex.end());//���մ�С��������
	auto it = std::unique(toDeleteTriIndex.begin(), toDeleteTriIndex.end());
	toDeleteTriIndex.erase(it, toDeleteTriIndex.end());//�Ƴ��ظ���
	for (size_t iTri = 0; iTri < triEleVec_.size(); iTri++)
	{
		auto target = triEleVec_.at(iTri).GetTriEleIndex();
		for (auto it : toDeleteTriIndex)
		{
			if (target == it)
			{
				triEleVec_.erase(triEleVec_.begin() + iTri);
				--iTri;
				break;
			}
		}
	}

	return;
}