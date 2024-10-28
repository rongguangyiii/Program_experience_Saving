#include "grid/include/node.h"

Node::Node(const Coord& c, const NodeId& id, NodeType type, NodeTag tag)
    : coord_(c), id_(id), type_(type), tag_(tag)
{
}

Node& Node::operator=(const Node& other) noexcept
{
    if (this != &other) {
        coord_ = other.coord_;
        id_ = other.id_;
        neighbors_ = other.neighbors_;
        type_ = other.type_;
        tag_ = other.tag_;
    }
    return *this;
}
