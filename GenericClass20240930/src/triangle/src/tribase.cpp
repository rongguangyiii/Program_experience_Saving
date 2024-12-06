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
#include "triangle/include/cgalTypes.h"
#include "triangle/include/triEle.h"
#include "triangle/include/triBase.h"
#include <memory>
#include <vector>

std::size_t EdgeHash::operator()(const std::shared_ptr<TriEdge>& k) const {
	std::size_t hash = 2166136261; // FNV offset basis  
	hash ^= k->vA;
	hash *= 16777619;
	hash ^= k->vB;
	hash *= 16777619;
	return hash;
}

void TriBase::addEdgeTableItem(const triElePtr_& triEle_ptr)
{
	for (const auto& edge : triEle_ptr->getTriEdgesRef())
	{
		if (edge->vA == edge->vB) {
			throw std::runtime_error("Caught an exception!");
		}
		edgeTable_[edge].push_back(triEle_ptr);
	}

}

TriBase::TriBase(std::vector<TriNode>& outterNode, std::vector<TriNode>& innerNode, const DataType& type)
	:dividType_(type), isoutFile(false)
{
	if (dividType_ == DataType::ring_ring){
		reorderNodes(outterNode);
		reorderNodes(innerNode);
	}
	else if (dividType_ == DataType::point_ring){
		reorderNodes(innerNode);
	}
	else if (dividType_ == DataType::point_line|| dividType_ == DataType::point_point){	
		//do nothing
	}
	else{
		throw std::runtime_error("Caught an exception!");
	}
	/*-----------------------------------------------------------------
	Generate TriNodeMap, each TriNode must be instantiated
	at the time of the point type, each point index must be unique;
	-----------------------------------------------------------------*/
	for (auto& per_it : outterNode)
	{
		per_it.setTriNodeType(TriNodeType::outter);
		const size_t& index = per_it.getindex();
		TriNodeMap_[index] = per_it;
		cgalPointsIndexVec_.emplace_back(index);
	}
	for (auto& per_it : innerNode)
	{
		per_it.setTriNodeType(TriNodeType::inner);
		const size_t& index = per_it.getindex();
		TriNodeMap_[index] = per_it;
		cgalPointsIndexVec_.emplace_back(index);
	}
}

void TriBase::insertCtrlPoly(std::vector<TriNode> nodes, TriNodeType nodetype)
{
	reorderNodes(nodes);
	for (auto& it : nodes)
	{
		it.setTriNodeType(nodetype);
		size_t index = it.getindex();
		TriNodeMap_[index] = it;
		cgalPointsIndexVec_.emplace_back(index);
	}
}
void TriBase::insertCtrlLine(std::vector<TriNode> nodes, TriNodeType nodetype)
{
	for (auto& it : nodes)
	{
		it.setTriNodeType(TriNodeType::inner);
		size_t index = it.getindex();
		TriNodeMap_[index] = it;
		cgalPointsIndexVec_.emplace_back(index);
	}
}
void TriBase::insertPoint(std::vector<TriNode> nodes, TriNodeType nodetype)
{
	for (auto& it : nodes)
	{
		it.setTriNodeType(nodetype);
		size_t index = it.getindex();
		TriNodeMap_[index] = it;
		cgalPointsIndexVec_.emplace_back(index);
	}
}

/*----------------------------------------------------------------------------------------------
	TO DEBUG
------------------------------------------------------------------------------------------------*/
TriBase::TriBase(const DataType& type)
	:dividType_(type), isoutFile(true)
{
	std::vector<std::vector<double>> poly_1_outter = { {1,2}, {3,2}, {7,2}, {7,4}, {8,7}, {5,8}, {2,6}, {3,4} };
	std::vector<std::vector<double>> poly_2_inner = { {4,4}, {6,3}, {6,6}, {4,6}, {5,5} };
	size_t index = 1;
	for (const auto& per_it : poly_1_outter)
	{
		TriNode Point_cdt{ per_it[0],per_it[1],0, index, TriNodeType::outter };
		TriNodeMap_[index] = Point_cdt;
		cgalPointsIndexVec_.emplace_back(index);
		index++;
	}
	for (const auto& per_it : poly_2_inner)
	{
		TriNode Point_cdt{ per_it[0],per_it[1],0, index, TriNodeType::inner };
		TriNodeMap_[index] = Point_cdt;
		cgalPointsIndexVec_.emplace_back(index);
		index++;
	}
}

void TriBase::genEdgeTable()
{
	edgeTable_.clear();
	for (auto& triEle : triEleVec_)
		addEdgeTableItem(triEle.second);

	return;
}

void TriBase::GenCompleteCell()
{
	for (auto& triEle : triEleVec_)
	{
		if (TriEleTag::used == triEle.second->getTriEleTag())
			continue;
		cell_.push_back(triEle.second->GetTriEleVertex());
	}

	return;
}

void TriBase::OutputConnectGrid2Tecplot(const std::vector<size_t> pointIndexVec, cellVec_& cellVec,
	const std::unordered_map<size_t, size_t>& vertex2tecplotIndex, const std::string title)const
{
	if (!isoutFile)return;
	std::string outfilename = "Data_" + title + "_grid" + ".dat";
	std::string basePath = "tempData";
	createFolder(basePath);

	std::ofstream tecplotFile(basePath + "/" + outfilename);
	tecplotFile << "TITLE = \"Triangle Output   LGYYYY\"\n";
	tecplotFile << "VARIABLES = \"X\", \"Y\"\n";

	tecplotFile << "ZONE T=\"Triangles\", N=" << pointIndexVec.size()
		<< ", E=" << cellVec.size() << ", F=FEPOINT, ET=QUADRILATERAL\n";
	for (const auto& vit : pointIndexVec)
	{
		const auto& p = index2Coor(vit);
		tecplotFile << p.x() << " " << p.y() << "\n";
	}

	for (const auto& cell : cellVec)
	{
		if (cell.size() == 3) {
			tecplotFile << vertex2tecplotIndex.at(cell[0]) + 1 << " " << vertex2tecplotIndex.at(cell[1]) + 1 << " "
				<< vertex2tecplotIndex.at(cell[2]) + 1 << " " << vertex2tecplotIndex.at(cell[0]) + 1 << std::endl;
		}
		else if (cell.size() == 4) { // �ı��ε�Ԫ
			tecplotFile << vertex2tecplotIndex.at(cell[0]) + 1 << " " << vertex2tecplotIndex.at(cell[1]) + 1 << " "
				<< vertex2tecplotIndex.at(cell[2]) + 1 << " " << vertex2tecplotIndex.at(cell[3]) + 1 << std::endl;
		}
		else {

		}
	}

	tecplotFile.close();

	return;
}

void TriBase::genCDT()
{
	std::vector<Point_2> cgal_outerPointsVec;
	std::vector<Point_2> cgal_innerPointsVec;
	std::vector<TriNode> poly_outerCoor;
	std::vector<TriNode> poly_innerCoor;

	for (const auto& perit : cgalPointsIndexVec_)
	{
		const auto& cur_trinode = TriNodeMap_.at(perit);
		if (cur_trinode.gettype() == TriNodeType::outter)
		{
			cgal_outerPointsVec.push_back({ cur_trinode.x() ,cur_trinode.y() });
			poly_outerCoor.emplace_back(cur_trinode);
		}
		else if (cur_trinode.gettype() == TriNodeType::inner)
		{
			cgal_innerPointsVec.push_back({ cur_trinode.x() ,cur_trinode.y() });
			poly_innerCoor.emplace_back(cur_trinode);
		}
		else
		{
			throw std::runtime_error("Caught an exception!");
		}
	}

	std::vector<Point_2> cgal_PointsVec(cgal_outerPointsVec.begin(), cgal_outerPointsVec.end());
	cgal_PointsVec.insert(cgal_PointsVec.end(), cgal_innerPointsVec.begin(), cgal_innerPointsVec.end());

	CDT cdt;
	std::unordered_map<size_t, size_t> vertex2tecplotIndex;
	std::map<Vertex_handle, size_t> vertex_to_index;
	for (size_t iPoint = 0; iPoint != cgalPointsIndexVec_.size(); ++iPoint)
	{
		Vertex_handle vh = cdt.insert(cgal_PointsVec.at(iPoint));
		vertex_to_index[vh] = cgalPointsIndexVec_.at(iPoint);
		vertex2tecplotIndex[cgalPointsIndexVec_.at(iPoint)] = iPoint;
	}
	//node idx
	vertex2tecplotIndex_.first = cgalPointsIndexVec_;
	//local idx
	vertex2tecplotIndex_.second = vertex2tecplotIndex;

	//indeart Constrain respectly
	if (dividType_ == DataType::ring_ring)
	{
		for (size_t i = 0; i < cgal_outerPointsVec.size(); ++i)
			cdt.insert_constraint(cgal_outerPointsVec[i], cgal_outerPointsVec[(i + 1) % cgal_outerPointsVec.size()]);
		for (size_t i = 0; i < cgal_innerPointsVec.size(); ++i)
			cdt.insert_constraint(cgal_innerPointsVec[i], cgal_innerPointsVec[(i + 1) % cgal_innerPointsVec.size()]);
	}
	else if (dividType_ == DataType::point_line|| dividType_ == DataType::point_ring)
	{
		for (size_t i = 0; i < cgal_innerPointsVec.size() - 1; ++i)
			cdt.insert_constraint(cgal_innerPointsVec[i], cgal_innerPointsVec[i + 1]);
	}
	else if (dividType_ == DataType::point_point)
	{ 
		//DO NOTHING
	}
	else
	{ 
		throw std::runtime_error("Data Type is wrong !");
	}

	size_t triIndex = 0;
	for (Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit)
	{
		size_t v0 = vertex_to_index[fit->vertex(0)];
		size_t v1 = vertex_to_index[fit->vertex(1)];
		size_t v2 = vertex_to_index[fit->vertex(2)];
		TriNode p0 = index2Coor(v0);
		TriNode p1 = index2Coor(v1);
		TriNode p2 = index2Coor(v2);
		//step 1 : 判断三角形的三个顶点是否在一条直线上,
		if (isCollinear(p0, p1, p2))
			continue;
		TriNode triCenter = (p0 += p1 += p2) /= 3.0;
		if (dividType_ == DataType::ring_ring) {
			if (isPointInPolygon_moreInner(triCenter, poly_innerCoor))continue;
			if (!isPointInPolygon_moreOuter(triCenter, poly_outerCoor))continue;
		}
		else if (dividType_ == DataType::point_ring) {
			if (isPointInPolygon_moreInner(triCenter, poly_innerCoor))continue;
			std::vector<TriNode> Nodes{ TriNodeMap_.at(v0) ,TriNodeMap_.at(v1),TriNodeMap_.at(v2) };
			if (judgeTriType(Nodes) == TriEleType::Outer_only)continue;
		}
		else if (dividType_ == DataType::point_line || dividType_ == DataType::point_point) {
			std::vector<TriNode> Nodes{ TriNodeMap_.at(v0) ,TriNodeMap_.at(v1),TriNodeMap_.at(v2) };
			if (judgeTriType(Nodes) == TriEleType::Outer_only)continue;
		}
		std::vector<TriNode> triNodes{ TriNodeMap_.at(v0) ,TriNodeMap_.at(v1),TriNodeMap_.at(v2) };
		TriEle triEle(triNodes, triIndex++);
		triEle.setTriCenter(triCenter);
		std::shared_ptr<TriEle> triEle_ptr = std::make_shared<TriEle>(triEle);
		addTriEle(triEle_ptr);
	}

	if (isoutFile) {
		std::vector< std::vector<size_t> > cellVec;
		for (const auto& perTri : triEleVec_)
			cellVec.push_back(perTri.second->GetTriEleVertex());
		OutputConnectGrid2Tecplot(vertex2tecplotIndex_.first, cellVec, vertex2tecplotIndex_.second, "Tri_ini");
	}
	return;
}

TriEleType TriBase::judgeTriType(std::vector<TriNode> triNodes) const
{
	size_t out_num = 0;
	size_t in_num = 0;
	for (auto& perit : triNodes) {
		if (perit.gettype() == TriNodeType::outter)
			out_num++;
		else if (perit.gettype() == TriNodeType::inner)
			in_num++;
		else
		{
			throw std::runtime_error("Caught an exception!");
		}
	}
	if (out_num == 3 && in_num == 0) {
		return TriEleType::Outer_only;
	}
	else if (out_num == 2 && in_num == 1) {
		return TriEleType::Outer;
	}
	else if (out_num == 1 && in_num == 2) {
		return TriEleType::Inner;
	}
	else if (out_num == 0 && in_num == 3) {
		return TriEleType::Inner_only;
	}
	else {
		throw std::runtime_error("Caught an exception!");
	}
}

void TriBase::filterTri()
{
	const auto& triEleVEcPtr = GetTriEleVec();
	for (const auto& perit : triEleVEcPtr)
	{
		auto& currentTri = perit.second;
		if (currentTri->getTriEleType() == TriEleType::Inner_only ||
			currentTri->getTriEleType() == TriEleType::Outer_only) {
			addNotIdealTri(*currentTri);
		}
	}
	return;
}

void TriBase::reconstructNotIdealTriGrid()
{
	filterTri();
	std::vector<TriEle>& notIdealTriVec = getNotIdealTriVec();
	if (0 == notIdealTriVec.size())
		return;

	genEdgeTable();
	std::vector<size_t> toDeleteTriIndex;
	for (size_t inot = 0; inot < notIdealTriVec.size(); inot++)
	{
		auto& currentTri = notIdealTriVec.at(inot);
		for (const auto& edge : currentTri.getTriEdgesRef())
		{
			auto it = getEdgeTable().find(edge);
			if (it->second.size() == 1)
				continue;

			const std::vector<std::shared_ptr<TriEle>>& triangles = it->second;
			for (auto& var : triangles)
				toDeleteTriIndex.push_back(var->GetTriEleIndex());
			reconstructNotIdealTriGrid(*edge, triangles);
		}
	}
	for (size_t cc : toDeleteTriIndex)
		removeTriEle(cc);

	if (isoutFile) {
		std::vector< std::vector<size_t> > cellVec;
		for (const auto& perTri : GetTriEleVec())
			cellVec.push_back(perTri.second->GetTriEleVertex());
		OutputConnectGrid2Tecplot(vertex2tecplotIndex_.first, cellVec, vertex2tecplotIndex_.second, "Tri_reconstruct");
	}
	return;
}

void TriBase::reconstructNotIdealTriGrid(const TriEdge& edge, const std::vector<triElePtr_>& tri_Ptr)
{
	std::vector<size_t> tri_1_vertex = tri_Ptr.at(0)->GetTriEleVertex();
	std::vector<size_t> tri_2_vertex = tri_Ptr.at(1)->GetTriEleVertex();
	size_t coindex_a = edge.vA;
	size_t coindex_b = edge.vB;
	size_t tri_1_diff = 0;
	size_t tri_2_diff = 0;
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
	std::vector<std::pair<TriNode, size_t >> newTri_1_pair;
	std::vector<std::pair<TriNode, size_t >> newTri_2_pair;
	TriNode average_1;
	TriNode average_2;
	for (auto it : newTri_1)
	{
		TriNode temcoor = index2Coor(it);
		newTri_1_pair.push_back(std::make_pair(temcoor, it));
		average_1 += (temcoor /= (double)newTri_1.size());
	}
	for (auto it : newTri_2)
	{
		TriNode temcoor = index2Coor(it);
		newTri_2_pair.push_back(std::make_pair(temcoor, it));
		average_2 += (temcoor /= (double)newTri_2.size());
	}
	reorderPointsIndex(newTri_1_pair, average_1);
	reorderPointsIndex(newTri_2_pair, average_2);
	newTri_1.clear();
	newTri_2.clear();
	for (auto& it : newTri_1_pair)
		newTri_1.push_back(it.second);
	for (auto& it : newTri_2_pair)
		newTri_2.push_back(it.second);

	std::vector<TriNode> triNodes_1;
	triNodes_1.emplace_back(TriNodeMap_.at(newTri_1.at(0)));
	triNodes_1.emplace_back(TriNodeMap_.at(newTri_1.at(1)));
	triNodes_1.emplace_back(TriNodeMap_.at(newTri_1.at(2)));
	TriEle triEle_1(triNodes_1, GetTriEleVec().size());
	triEle_1.setTriCenter(average_1);
	const auto& current_triptr = std::make_shared<TriEle>(triEle_1);
	addTriEle(current_triptr);
	addEdgeTableItem(current_triptr);

	std::vector<TriNode> triNodes_2;
	triNodes_2.emplace_back(TriNodeMap_.at(newTri_2.at(0)));
	triNodes_2.emplace_back(TriNodeMap_.at(newTri_2.at(1)));
	triNodes_2.emplace_back(TriNodeMap_.at(newTri_2.at(2)));
	TriEle triEle_2(triNodes_2, GetTriEleVec().size());
	triEle_2.setTriCenter(average_2);
	const auto& another_ptr = std::make_shared<TriEle>(triEle_2);
	addTriEle(another_ptr);
	addEdgeTableItem(another_ptr);

	return;
}

void TriBase::genTriAndQuad()
{
	if (0 == edgeTable_.size())
		genEdgeTable();

	for (const auto& perit : edgeTable_)
	{
		auto& cur_edge = perit.first;
		auto& cur_triVec = perit.second;
		if (cur_triVec.size() == 0 || cur_triVec.size() >= 3) {
			throw std::runtime_error("Caught an exception!");
		}
		if (cur_triVec.size() == 1)
			continue;
		if (cur_triVec[0]->getTriEleTag() == TriEleTag::used ||
			cur_triVec[1]->getTriEleTag() == TriEleTag::used)
			continue;
		if (cur_triVec[0]->getTriEleType() == cur_triVec[1]->getTriEleType())
			continue;

		genTriAndQuad(*cur_edge, cur_triVec);

	}//end for edge

	GenCompleteCell();
	if (isoutFile){
		OutputConnectGrid2Tecplot(vertex2tecplotIndex_.first, getCell(), vertex2tecplotIndex_.second, "Tri_quad_final");
	}

	return;
}

void TriBase::genTriAndQuad(const TriEdge& edge, const std::vector<triElePtr_>& triangles)
{
	std::vector<size_t> tri_1_vertex = triangles.at(0)->GetTriEleVertex();
	std::vector<size_t> tri_2_vertex = triangles.at(1)->GetTriEleVertex();
	size_t coindex_a = edge.vA;
	size_t coindex_b = edge.vB;
	size_t tri_1_diff = 0;
	size_t tri_2_diff = 0;
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
	std::vector<std::pair<TriNode, size_t >> quad_pair;
	TriNode average;
	for (auto it : quad)
	{
		TriNode temcoor = index2Coor(it);
		quad_pair.push_back(std::make_pair(temcoor, it));
		average += (temcoor /= (double)quad.size());
	}

	reorderPointsIndex(quad_pair, average);
	quad.clear();
	for (auto& it : quad_pair)
		quad.push_back(it.second);

	// The first time getCell is cell_ is empty, add cell_ here
	getCell().push_back(quad);

	//TriEleTag::used
	for (auto& var : triangles)
	{
		if (var->getTriEleTag() == TriEleTag::used) {
			throw std::runtime_error("Caught an exception!");
		}
		var->setTriEleTag(TriEleTag::used);
	}

	return;
}

void TriBase::gencell(const triEleMap_& content)
{
	for (auto& triEle : content)
		cell_.push_back(triEle.second->GetTriEleVertex());
	return;
}

void TriBase::boostMesh()
{
	genCDT();
	//reconstructNotIdealTriGrid();//No need!
	std::string method = "CDT";
	//method = "HYB";//DEBUG
	if (method == "CDT")
		gencell(GetTriEleVec());
	else if (method == "HYB")
		genTriAndQuad();
	else
		throw std::runtime_error("TriEle not found in triEleVec_");
	isRemarkBack();
	return;
}

void TriBase::isRemarkBack()
{
	for (auto& per : TriNodeMap_) {
		auto& idx = per.first;
		bool ishave = false;
		for (auto& perEle : cell_) {
			for (auto& it : perEle) {
				if (it == idx) {
					ishave = true;
					break;
				}
			}
			if (ishave)break;
		}
		if (!ishave) {
			reMarkIdx_.push_back(idx);
		}
	}
}

void TriBase::addTriEle(const triElePtr_& triEle)
{
	size_t index = triEle->GetTriEleIndex();
	triEleVec_.insert_or_assign(index, triEle);
};

void TriBase::removeTriFromEdgeTable(const triElePtr_& triEle_ptr) {
	for (const auto& edge : triEle_ptr->getTriEdgesRef()) {
		auto it = edgeTable_.find(edge);
		if (it != edgeTable_.end()) {
			auto& triEleList = it->second;
			triEleList.erase(std::remove(triEleList.begin(), triEleList.end(), triEle_ptr), triEleList.end());
			if (triEleList.empty()) {
				edgeTable_.erase(it);
			}
		}
	}
}

void TriBase::removeTriEle(size_t triEleKey) {
	auto it = triEleVec_.find(triEleKey);
	if (it != triEleVec_.end()) {
		auto& triEle_ptr = it->second;
		removeTriFromEdgeTable(triEle_ptr);
		triEleVec_.erase(it);
	}
	else {
		throw std::runtime_error("TriEle not found in triEleVec_");
	}
}

//--------------------------------------------------------------------------------------------------------------------------
void TriBase::reorderPointsIndex(std::vector<std::pair<TriNode, size_t>>& PointLists, const TriNode& modelCenter)
{
	std::size_t wallNodeNum = PointLists.size();
	std::vector<std::pair<double, size_t>> angleAndIndex;
	angleAndIndex.resize(wallNodeNum);
	for (size_t iNode = 0; iNode < wallNodeNum; ++iNode)
	{
		auto& perNode = PointLists.at(iNode).first;
		TriNode vec = perNode - modelCenter;
		double angle = atan2(vec.y(), vec.x());
		angleAndIndex[iNode] = std::make_pair(angle, iNode);
	}
	std::sort(angleAndIndex.begin(), angleAndIndex.end(),
		[](const std::pair<double, size_t>& a, const std::pair<double, size_t>& b) {return a.first < b.first; });

	std::vector<size_t> reorderedIndices;
	std::vector<TriNode> NodeVec_temp;
	std::vector<std::pair<TriNode, size_t>> NodePairVec_temp2;
	NodeVec_temp.resize(wallNodeNum);
	reorderedIndices.resize(wallNodeNum);
	for (size_t i = 0; i < wallNodeNum; ++i)
	{
		NodeVec_temp[i] = PointLists[angleAndIndex[i].second].first;
		reorderedIndices[i] = PointLists[angleAndIndex[i].second].second;
		NodePairVec_temp2.push_back(std::make_pair(NodeVec_temp[i], reorderedIndices[i]));
	}
	PointLists = NodePairVec_temp2;
}

bool TriBase::isCollinear(const TriNode& p1, const TriNode& p2, const TriNode& p3)
{
	constexpr double epsilon = std::numeric_limits<double>::epsilon();
	double dx1 = p2.x() - p1.x();
	double dy1 = p2.y() - p1.y();
	double dx2 = p3.x() - p1.x();
	double dy2 = p3.y() - p1.y();

	if (std::abs(dx1 * dy2 - dx2 * dy1) < epsilon) {
		return true;
	}
	if (std::abs(dx1) < epsilon && std::abs(dx2) < epsilon) {
		return std::abs(dy1) < epsilon && std::abs(dy2) < epsilon;
	}

	return false;
}

void TriBase::createFolder(const std::string& basePath)const
{
	std::filesystem::path folderPath = basePath;

	if (!std::filesystem::exists(folderPath)) {
		std::filesystem::create_directories(folderPath);
	}
}

bool TriBase::isPointInPolygon_moreInner(const TriNode& p, const std::vector<TriNode>& polygon) {
	const double EPSILON = 1e-10;
	bool inside = false;
	size_t j = polygon.size() - 1;

	for (size_t i = 0; i < polygon.size(); i++) {
		// Check if point is on vertex
		if ((std::fabs(p.x() - polygon[i].x()) < EPSILON) && (std::fabs(p.y() - polygon[i].y()) < EPSILON)) {
			return true; // Point is on a vertex
		}

		// Check if point is on edge
		double minX = std::min(polygon[i].x(), polygon[j].x());
		double maxX = std::max(polygon[i].x(), polygon[j].x());
		double minY = std::min(polygon[i].y(), polygon[j].y());
		double maxY = std::max(polygon[i].y(), polygon[j].y());

		if ((minX <= p.x() && p.x() <= maxX) && (minY <= p.y() && p.y() <= maxY)) {
			double dx = polygon[j].x() - polygon[i].x();
			double dy = polygon[j].y() - polygon[i].y();
			if (std::fabs(dy) < EPSILON) { // horizontal line
				if (std::fabs(p.y() - polygon[i].y()) < EPSILON) {
					return true; // Point is on horizontal edge
				}
			}
			else if (std::fabs(dx) < EPSILON) { // vertical line
				if (std::fabs(p.x() - polygon[i].x()) < EPSILON) {
					return true; // Point is on vertical edge
				}
			}
			else {
				double slope = dy / dx;
				double intercept = polygon[i].y() - slope * polygon[i].x();
				if (std::fabs(p.y() - (slope * p.x() + intercept)) < EPSILON) {
					return true; // Point is on non-vertical edge
				}
			}
		}

		// Check if point is inside using ray-casting algorithm
		if (((polygon[i].y() <= p.y()) && (p.y() < polygon[j].y())) ||
			((polygon[j].y() <= p.y()) && (p.y() < polygon[i].y()))) {
			double intersectionX = polygon[i].x() + (p.y() - polygon[i].y()) * (polygon[j].x() - polygon[i].x()) / (polygon[j].y() - polygon[i].y());
			if (intersectionX < p.x()) {
				inside = !inside;
			}
		}

		j = i;
	}
	return inside;
}


bool TriBase::isPointInPolygon_moreOuter(const TriNode& p, const std::vector<TriNode>& polygon) {
	const double EPSILON = 1e-10;
	bool inside = false;
	size_t j = polygon.size() - 1;

	for (size_t i = 0; i < polygon.size(); i++) {
		// Check if point is on vertex
		if ((std::fabs(p.x() - polygon[i].x()) < EPSILON) && (std::fabs(p.y() - polygon[i].y()) < EPSILON)) {
			return false; // Point is on a vertex
		}

		// Check if point is on edge
		double minX = std::min(polygon[i].x(), polygon[j].x());
		double maxX = std::max(polygon[i].x(), polygon[j].x());
		double minY = std::min(polygon[i].y(), polygon[j].y());
		double maxY = std::max(polygon[i].y(), polygon[j].y());

		if ((minX <= p.x() && p.x() <= maxX) && (minY <= p.y() && p.y() <= maxY)) {
			double dx = polygon[j].x() - polygon[i].x();
			double dy = polygon[j].y() - polygon[i].y();
			if (std::fabs(dy) < EPSILON) { // horizontal line
				if (std::fabs(p.y() - polygon[i].y()) < EPSILON) {
					return false; // Point is on horizontal edge
				}
			}
			else if (std::fabs(dx) < EPSILON) { // vertical line
				if (std::fabs(p.x() - polygon[i].x()) < EPSILON) {
					return false; // Point is on vertical edge
				}
			}
			else {
				double slope = dy / dx;
				double intercept = polygon[i].y() - slope * polygon[i].x();
				if (std::fabs(p.y() - (slope * p.x() + intercept)) < EPSILON) {
					return false; // Point is on non-vertical edge
				}
			}
		}

		// Check if point is inside using ray-casting algorithm
		if (((polygon[i].y() <= p.y()) && (p.y() < polygon[j].y())) ||
			((polygon[j].y() <= p.y()) && (p.y() < polygon[i].y()))) {
			double intersectionX = polygon[i].x() + (p.y() - polygon[i].y()) * (polygon[j].x() - polygon[i].x()) / (polygon[j].y() - polygon[i].y());
			if (intersectionX < p.x()) {
				inside = !inside;
			}
		}

		j = i;
	}
	return inside;
}

void TriBase::reorderNodes(std::vector<TriNode>& elements)
{

	std::size_t elementNum = elements.size();
	TriNode modelCenter{ 0,0,0 };
	for (size_t iPoint = 0; iPoint < elementNum; ++iPoint)
	{
		TriNode per = elements[iPoint];
		modelCenter += per;
	}
	modelCenter /= (double)elementNum;


	std::vector<std::pair<double, size_t>> angleAndIndex;
	angleAndIndex.reserve(elementNum);
	for (size_t i = 0; i < elementNum; ++i)
	{
		auto& element = elements[i];
		TriNode vec = element - modelCenter;
		double angle = atan2(vec.y(), vec.x());
		angleAndIndex.emplace_back(angle, i);
	}

	std::sort(angleAndIndex.begin(), angleAndIndex.end(),
		[](const std::pair<double, size_t>& a, const std::pair<double, size_t>& b) { return a.first < b.first; });

	std::vector<TriNode> elementsTemp(elementNum);
	for (size_t i = 0; i < elementNum; ++i)
	{
		elementsTemp[i] = elements[angleAndIndex[i].second];
	}

	elements = std::move(elementsTemp);
}