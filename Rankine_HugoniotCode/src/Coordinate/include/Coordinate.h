/*---------------------------------------------------------------------------
	ZaRan	-	A Totally Automatic CFD Software
	Copyright (C) ,Since 2020
-------------------------------------------------------------------------------
License
	This file is part of ZaRan.

!	@file		gird.h
!	@brief	the purpose of this file.
!	@author	Chen Jie.
\*---------------------------------------------------------------------------*/
#pragma once
#include<vector>
//#include"Utility/include/math.h"
using std::vector;
/*
几何基础名称空间
由点、线、面基础类构成
*/

enum class CoorType
{
	ordinary,
	patch,
};

class Coordinate
{
public:
	Coordinate(double x = 0, double y = 0, double z = 0, CoorType type = CoorType::ordinary) :x_(x), y_(y), z_(z), coorType_(type) {}
	Coordinate(const Coordinate& pt) :x_(pt.x_), y_(pt.y_), z_(pt.z_), coorType_(pt.coorType_) {}
	const double& x() const { return x_; }
	const double& y() const { return y_; }
	const double& z() const { return z_; }
	void SetCoord(const double& x = 0, const double &y = 0, const double &z = 0);
	bool operator ==(const Coordinate& pt)const
	{
		double delta = 1e-14;
		return fabs(x_ - pt.x_) < delta && fabs(y_ - pt.y_) < delta && fabs(z_ - pt.z_) < delta;
	};
	double norm() const;
	void setCoorType(CoorType type) { coorType_ = type; }
	CoorType getCoorType() const { return coorType_; }
	Coordinate& operator+=(const Coordinate& pt);
	Coordinate& operator-=(const Coordinate& pt);
	Coordinate& operator/=(double scalar);
	Coordinate operator-(const Coordinate& pt)const;
	Coordinate operator+(const Coordinate& pt)const;
	bool operator<(const Coordinate& other) const;
private:
	double x_, y_, z_;
	CoorType coorType_;
};

struct CoorHash {
	std::size_t operator()(const Coordinate& pt) const;
};

double distance(const Coordinate& pt1, const Coordinate& pt2);
double distance(double x1, double y1, double z1, double x2, double y2, double z2);