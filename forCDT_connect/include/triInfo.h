#pragma once
#include <vector> 
#include <unordered_map>

struct Point {
	double x, y;

	Point(double x = 0, double y = 0) : x(x), y(y) {}

	operator const Point& () const {
		return *this;
	}
	//重载-运算符
	Point operator-(const Point& p) const {
		return Point(x - p.x, y - p.y);
	}
};


enum class EdgeType
{
	inner,
	boundary,
	unset,
};
// 定义边结构，使用顶点索引表示
class TriEdge
{
private:

public:
	size_t v1, v2;
	EdgeType edgeType_;
public:
	TriEdge();
	TriEdge(int v1, int v2);

	// 为了使用Edge作为unordered_map的键，需要重载==和hash函数
	bool operator==(const TriEdge& other) const {
		return v1 == other.v1 && v2 == other.v2;
	}
	void setEdgeType(EdgeType edgeType) { edgeType_ = edgeType; }
	EdgeType getEdgeType() const { return edgeType_; }

};

struct EdgeHash {
	std::size_t operator()(const TriEdge& k) const;
};

class TriEle
{
private:
	size_t EleIndex_;
	size_t a_, b_, c_;
	Point baryCenter_;

public:
	TriEle(size_t a, size_t b, size_t c, size_t index);
	TriEle();
	std::vector< size_t> GetTriEleVertex() const;
	size_t GetTriEleIndex() const { return EleIndex_; };
	void reSetTriEleIndex(size_t currentTriIndex) { EleIndex_ = currentTriIndex;};
	Point getBaryCenter() const { return baryCenter_; }
	void setBaryCenter(const Point& baryCenter) { baryCenter_ = baryCenter; }
};

class TriBase
{
	std::vector<TriEle> triEleVec_;
	std::vector<TriEle> notIdealTriVec_;
	std::vector<size_t> tecplotIndex_;
	std::vector<Point> all_points_;
	std::vector<Point> poly_outer;
	std::vector<Point> poly_inner;
public:

	TriBase() {};
	TriEle GetTriEle(size_t index) const;
	std::vector<TriEle>& GetTriEleVec() { return triEleVec_; };
	void reSetAllTriEleIndex();
	void AddTriEle(TriEle triEle);
	int genCDT();
	void filterTriGrid();
	void checkTriGrid();
	void OutputTriGrid2Tecplot(const std::vector<size_t> pointIndexVec, std::vector<TriEle> triVec, const std::string title);
	void addTri2EdgeTable(TriEle& triEle);
	std::vector<TriEle> getNotIdealTriVec() const { return notIdealTriVec_; }
	void addNotIdealTri(TriEle triEle);
	void reconstructNotIdealTriGrid();
	void reconstructNotIdealTriGrid(const TriEdge& edge, const std::vector<TriEle>& triangles);
	void DeleteTri(std::vector<size_t> toDeleteTriIndex);
	std::vector<size_t> returnTecplotIndex();

public:
	// 边表使用unordered_map来存储，键是Edge，值是包含该边的三角形列表
	using EdgeTable = std::unordered_map<TriEdge, std::vector<TriEle>, EdgeHash>;
	EdgeTable edgeTable_;
};

bool isPointInPolygon(const Point& p, const std::vector<Point>& polygon);
void reorderPointsIndex(std::vector<std::pair<Point, size_t>>& PointLists, const Point& modelCenter);


