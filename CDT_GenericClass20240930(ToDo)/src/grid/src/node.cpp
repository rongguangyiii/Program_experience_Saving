#include "grid/include/node.h"

Node::Node(const Coord& c, const NodeId& id, NodeType typ, NodeTag ta)
	:coord_(c), id_(id), type_(typ), tag_(ta)
{

}
//写一个重载=的函数
