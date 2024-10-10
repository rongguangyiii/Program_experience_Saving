#pragma once
#include "grid/include/Coord.h"
#include <vector>

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
	Coord coord_;//����
	NodeId id_;  //���
	std::vector<NodeId> neighbors_;//�ھ�,Ŀǰÿ������ھӻ�û�д���
	NodeType type_;
	NodeTag tag_;
public:
	//����type_
	NodeType& getType() { return type_; }
	NodeTag& getTag() { return tag_; }
	const NodeId& getid()const { return id_; }
	const Coord& getCoord()const { return coord_; }

};
using NodePtr = std::shared_ptr<Node>;
