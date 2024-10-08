#pragma once
#include <vector>
#include <iomanip>
#include <unordered_map>
#include <set>

enum class PointType
{
	unset,
	inner,
	b_notcal,
	b_Lbound,
	b_Rbound,
	b_wallbound,
	b_filedbound,
	s_filedbound
};

// 结构体定义来唯一标识一个点
struct PointId
{
	size_t gridId; // 网格编号
	size_t nodeId; // 节点编号
	PointId(size_t grid_id, size_t node_id) : gridId(grid_id), nodeId(node_id) {}
	bool operator==(const PointId& other) const
	{
		return gridId == other.gridId && nodeId == other.nodeId;
	}
};

// 定义自定义哈希函数以支持 std::unordered_map
struct PointIdHash
{
	std::size_t operator()(const PointId& id) const
	{
		return std::hash<size_t>()(id.gridId) ^ (std::hash<size_t>()(id.nodeId) << 1);
	}
};
//只针对背景网格点用于查找吐出的背景网格点
enum class PointTag
{
	surface,
	shade,
	spit
};

struct Coord
{
	double x, y, z;
	PointType type_;
	PointId id_;
	PointTag tag_;
	//Coord GetCoord()
	// z 默认值为 0.0，使得结构体可以用于二维和三维点
	Coord(PointId id, double x, double y, double z = 0.0) : x(x), y(y), z(z), id_(id), type_(PointType::unset), tag_(PointTag::surface) {}
	Coord(double x, double y, double z = 0.0) : x(x), y(y), z(z), type_(PointType::unset), id_(0, 0), tag_(PointTag::surface) {}
	Coord() : x(0), y(0), z(0), type_(PointType::unset), id_(0, 0), tag_(PointTag::surface) {}

	bool operator<(const Coord& other) const
	{
		return std::tie(x, y) < std::tie(other.x, other.y); // 只比较 x 和 y
	}

	bool operator==(const Coord& other) const
	{
		return std::fabs(x - other.x) < 1e-9 &&
			std::fabs(y - other.y) < 1e-9 &&
			std::fabs(z - other.z) < 1e-9;
	}

	void operator=(const Coord& other)
	{
		this->x = other.x;
		this->y = other.y;
		this->z = other.z;
		this->type_ = other.type_;
	}
	// 重载加法运算符
	Coord& operator+=(const Coord& other)
	{
		x += other.x;
		y += other.y;
		z += other.z;
		return *this;
	}

	// 重载除法运算符
	Coord& operator/=(double scalar)
	{
		if (scalar != 0.0) {
			x /= scalar;
			y /= scalar;
			z /= scalar;
		}
		return *this;
	}
	// 重载减法运算符
	Coord& operator-=(const Coord& other)
	{
		x -= other.x;
		y -= other.y;
		z -= other.z;
		return *this;
	}

};

using CoordPtr = std::shared_ptr<Coord>;

enum class Elementype
{
	out,
	Notout
};

class Element
{
public:
	Element(const std::vector<size_t>& idvec) : connectivty_(idvec), eleType_(Elementype::out) {}
	Element(size_t n1, size_t n2, size_t n3, size_t n4)
	{
		connectivty_.push_back(n1);
		connectivty_.push_back(n2);
		connectivty_.push_back(n3);
		connectivty_.push_back(n4);
		eleType_ = Elementype::out;
	}
	Element(size_t n1, size_t n2, size_t n3)
	{
		connectivty_.push_back(n1);
		connectivty_.push_back(n2);
		connectivty_.push_back(n3);
		eleType_ = Elementype::out;
	}
	~Element() {}
	const std::vector<size_t>& getConnectivty() const { return connectivty_; }
	void setConnectivty(std::vector<size_t> idvec) { connectivty_ = idvec; }
	Elementype getEleType() const { return eleType_; }
	void setEleType(Elementype type) { eleType_ = type; }

private:
	std::vector<size_t> connectivty_;
	Elementype eleType_;
};
