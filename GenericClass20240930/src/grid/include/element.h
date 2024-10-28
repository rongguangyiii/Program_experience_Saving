#pragma once
#include <vector>

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