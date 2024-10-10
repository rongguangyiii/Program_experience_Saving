#pragma once
#include "grid/include/Coord.h"
#include <vector>

struct NodeId
{
	size_t gridId_;  // 网格编号
	size_t pointId_; // 节点编号
	NodeId(size_t grid_id, size_t point_id) : gridId_(grid_id), pointId_(point_id) {}
	NodeId() : gridId_(0), pointId_(0) {}
	bool operator==(const NodeId& other) const
	{
		return gridId_ == other.gridId_ && pointId_ == other.pointId_;
	}
};

// 定义自定义哈希函数以支持 std::unordered_map
struct NodeIdHash
{
	std::size_t operator()(const NodeId& id) const
	{
		return std::hash<size_t>()(id.gridId_) ^ (std::hash<size_t>()(id.pointId_) << 1);
	}
};
enum class NodeType
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

enum class NodeTag
{
	surface,
	shade,
	spit
};

class Node
{
public:
	Node() :coord_(Coord(0, 0, 0)), id_(NodeId(0, 0)), type_(NodeType::unset), tag_(NodeTag::surface) {};
	Node(const Coord& c, const NodeId& id, NodeType typ = NodeType::unset, NodeTag ta = NodeTag::surface);
	~Node() {};

private:
	Coord coord_;//坐标
	NodeId id_;  //编号
	std::vector<NodeId> neighbors_;//邻居,目前每个点的邻居还没有处理
	NodeType type_;
	NodeTag tag_;
public:
	//返回type_
	NodeType& getType() { return type_; }
	NodeTag& getTag() { return tag_; }
	const NodeId& getid()const { return id_; }
	const Coord& getCoord()const { return coord_; }

};
using NodePtr = std::shared_ptr<Node>;
