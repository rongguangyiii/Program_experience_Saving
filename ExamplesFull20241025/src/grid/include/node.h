#pragma once
#include "grid/include/Coord.h"
#include <vector>
#include <memory>

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

// 使用 enum class 替代 enum
enum class NodeType {
    unset,
    inner,
    BNotCal,
    BLBound,
    BRBound,
    BWallBound,
    b_filedbound,
    s_filedbound
};

enum class NodeTag {
    Surface,
    Shade,
    Spit
};

class Node {
public:
    // 使用默认构造函数
    Node() = default;
    // 使用初始化列表和默认参数
    Node(const Coord& c, const NodeId& id, NodeType type = NodeType::unset, NodeTag tag = NodeTag::Surface);
    
    // 使用默认析构函数
    ~Node() = default;

    // 使用 const 引用返回只读数据
    NodeType& getType()  { return type_; }
    const NodeTag& getTag() const { return tag_; }
    const NodeId& getid() const { return id_; }
    const Coord& getCoord() const { return coord_; }
    const std::vector<NodeId>& getNeighbors() const { return neighbors_; }

    // 提供修改数据的方法
    void setType(NodeType type) { type_ = type; }
    void setTag(NodeTag tag) { tag_ = tag; }
    void addNeighbor(const NodeId& neighbor) { neighbors_.push_back(neighbor); }

    // 使用 noexcept 标记不抛出异常的函数
    Node& operator=(const Node& other) noexcept;

private:
    Coord coord_;
    NodeId id_;
    std::vector<NodeId> neighbors_;
    NodeType type_ = NodeType::unset;
    NodeTag tag_ = NodeTag::Surface;
};

// 使用 std::shared_ptr
using NodePtr = std::shared_ptr<Node>;
