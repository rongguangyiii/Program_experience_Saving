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

#include "triangle/include/trinode.h"
#include <stdexcept>
#include <cmath>

void TriNode::SetCoord(const double& x, const double& y, const double& z)
{
	x_ = x;
	y_ = y;
	z_ = z;
}
bool TriNode::operator ==(const TriNode& pt)const
{
	double delta = 1e-14;
	return fabs(x_ - pt.x_) < delta && fabs(y_ - pt.y_) < delta && fabs(z_ - pt.z_) < delta;
};
bool TriNode::operator<(const TriNode& other) const {
	if (x_ != other.x_) return x_ < other.x_;
	if (y_ != other.y_) return y_ < other.y_;
	return z_ < other.z_;
}
TriNode& TriNode::operator+=(const TriNode& pt)
{
	x_ += pt.x();
	y_ += pt.y();
	z_ += pt.z();
	return *this;
}
TriNode& TriNode::operator-=(const TriNode& pt)
{
	x_ -= pt.x();
	y_ -= pt.y();
	z_ -= pt.z();
	return *this;
}
TriNode TriNode::operator-(const TriNode& pt) const
{
	return TriNode(x_ - pt.x(), y_ - pt.y(), z_ - pt.z());
}
TriNode TriNode::operator+(const TriNode& pt) const
{
	return TriNode(x_ + pt.x(), y_ + pt.y(), z_ + pt.z());
}
TriNode& TriNode::operator/=(double scalar)
{
	if (scalar != 0)
	{
		x_ /= scalar;
		y_ /= scalar;
		z_ /= scalar;
	}
	else
	{
		// 处理除数为零的情况  
		throw std::invalid_argument("Cannot divide by zero.");
	}
	return *this;
}
double TriNode::norm() const
{
	return std::sqrt(x_ * x_ + y_ * y_);
}
