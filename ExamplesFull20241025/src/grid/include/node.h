#pragma once
#include "grid/include/Coord.h"
#include <vector>
#include <memory>

struct NodeId
{
	size_t gridId_;  // ������
	size_t pointId_; // �ڵ���
	NodeId(size_t grid_id, size_t point_id) : gridId_(grid_id), pointId_(point_id) {}
	NodeId() : gridId_(0), pointId_(0) {}
	bool operator==(const NodeId& other) const
	{
		return gridId_ == other.gridId_ && pointId_ == other.pointId_;
	}
};

// �����Զ����ϣ������֧�� std::unordered_map
struct NodeIdHash
{
	std::size_t operator()(const NodeId& id) const
	{
		return std::hash<size_t>()(id.gridId_) ^ (std::hash<size_t>()(id.pointId_) << 1);
	}
};

// ʹ�� enum class ��� enum
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
    // ʹ��Ĭ�Ϲ��캯��
    Node() = default;
    // ʹ�ó�ʼ���б��Ĭ�ϲ���
    Node(const Coord& c, const NodeId& id, NodeType type = NodeType::unset, NodeTag tag = NodeTag::Surface);
    
    // ʹ��Ĭ����������
    ~Node() = default;

    // ʹ�� const ���÷���ֻ������
    NodeType& getType()  { return type_; }
    const NodeTag& getTag() const { return tag_; }
    const NodeId& getid() const { return id_; }
    const Coord& getCoord() const { return coord_; }
    const std::vector<NodeId>& getNeighbors() const { return neighbors_; }

    // �ṩ�޸����ݵķ���
    void setType(NodeType type) { type_ = type; }
    void setTag(NodeTag tag) { tag_ = tag; }
    void addNeighbor(const NodeId& neighbor) { neighbors_.push_back(neighbor); }

    // ʹ�� noexcept ��ǲ��׳��쳣�ĺ���
    Node& operator=(const Node& other) noexcept;

private:
    Coord coord_;
    NodeId id_;
    std::vector<NodeId> neighbors_;
    NodeType type_ = NodeType::unset;
    NodeTag tag_ = NodeTag::Surface;
};

// ʹ�� std::shared_ptr
using NodePtr = std::shared_ptr<Node>;
