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
#include"triangle/include/triCoord.h"
#include<algorithm>
#include <cmath>
#include <stdexcept>

void TriCoord::SetCoord(const double& x, const double& y, const double& z)
{
	x_ = x;
	y_ = y;
	z_ = z;
}

bool TriCoord::operator<(const TriCoord& other) const {
	if (x_ != other.x_) return x_ < other.x_;
	if (y_ != other.y_) return y_ < other.y_;
	return z_ < other.z_;
}

TriCoord& TriCoord::operator+=(const TriCoord& pt)
{
	x_ += pt.x();
	y_ += pt.y();
	z_ += pt.z();
	return *this;
}
TriCoord& TriCoord::operator-=(const TriCoord& pt)
{
	x_ -= pt.x();
	y_ -= pt.y();
	z_ -= pt.z();
	return *this;
}
TriCoord TriCoord::operator-(const TriCoord& pt) const
{
	return TriCoord(x_ - pt.x(), y_ - pt.y(), z_ - pt.z());
}
TriCoord TriCoord::operator+(const TriCoord& pt) const
{
	return TriCoord(x_ + pt.x(), y_ + pt.y(), z_ + pt.z());
}
TriCoord& TriCoord::operator/=(double scalar)
{
	if (scalar != 0)
	{
		x_ /= scalar;
		y_ /= scalar;
		z_ /= scalar;
	}
	else
	{
		// �������Ϊ������  
		throw std::invalid_argument("Cannot divide by zero.");
	}
	return *this;
}

double TriCoord::norm() const
{
	return std::sqrt(x_ * x_ + y_ * y_);
}

size_t CoorHash::operator()(const TriCoord& pt) const {

	// ֱ��ʹ�ø����������й�ϣ������ܻᵼ�²�ͬ�ĸ�����������ͬ�Ĺ�ϣֵ����������ĸ��������н�Ȼ��ͬ�Ĺ�ϣֵ�� 
	// Ϊ�������ڸ��������ȴ�����Ӱ�죬ʹ�ù̶��ľ�������������ֵ 
	const double precision = 1e-7; // ���磬ʹ��С�����8λ�ľ��ȣ������ļ�����ֻ��С�����7λ.  
	double qx = std::round(pt.x() / precision) * precision;
	double qy = std::round(pt.y() / precision) * precision;
	double qz = std::round(pt.z() / precision) * precision;

	// ʹ�������������ֵ�������ϣ  
	auto hash1 = std::hash<double>{}(qx);
	auto hash2 = std::hash<double>{}(qy);
	auto hash3 = std::hash<double>{}(qz);

	// ���������ϣֵ�Բ������յĹ�ϣ  
	return hash1 ^ (hash2 << 1) ^ (hash3 << 2);
}
