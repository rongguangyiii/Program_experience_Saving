#pragma once
#include "grid/include/Coord.h"

bool Coord::operator<(const Coord& other) const
{
	return std::tie(x_, y_) < std::tie(other.x_, other.y_); // 只比较 x 和 y
}

bool Coord::operator==(const Coord& other) const
{
	return std::fabs(x_ - other.x_) < 1e-9 &&
		std::fabs(y_ - other.y_) < 1e-9 &&
		std::fabs(z_ - other.z_) < 1e-9;
}

void Coord::operator=(const Coord& other)
{
	this->x_ = other.x_;
	this->y_ = other.y_;
	this->z_ = other.z_;
}
// 重载加法运算符
Coord& Coord::operator+=(const Coord& other)
{
	x_ += other.x_;
	y_ += other.y_;
	z_ += other.z_;
	return *this;
}

// 重载除法运算符
Coord& Coord::operator/=(double scalar)
{
	if (scalar != 0.0) {
		x_ /= scalar;
		y_ /= scalar;
		z_ /= scalar;
	}
	return *this;
}
// 重载减法运算符
Coord& Coord::operator-=(const Coord& other)
{
	x_ -= other.x_;
	y_ -= other.y_;
	z_ -= other.z_;
	return *this;
}
