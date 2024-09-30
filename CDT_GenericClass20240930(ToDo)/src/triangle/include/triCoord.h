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
#pragma once
#include<vector>
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

class TriCoord
{
public:
	TriCoord(double x = 0, double y = 0, double z = 0, CoorType type = CoorType::ordinary) :x_(x), y_(y), z_(z), coorType_(type) {}
	TriCoord(const TriCoord& pt) :x_(pt.x_), y_(pt.y_), z_(pt.z_), coorType_(pt.coorType_) {}
	const double& x() const { return x_; }
	const double& y() const { return y_; }
	const double& z() const { return z_; }
	void SetCoord(const double& x = 0, const double &y = 0, const double &z = 0);
	bool operator ==(const TriCoord& pt)const
	{
		double delta = 1e-14;
		return fabs(x_ - pt.x_) < delta && fabs(y_ - pt.y_) < delta && fabs(z_ - pt.z_) < delta;
	};
	double norm() const;
	void setCoorType(CoorType type) { coorType_ = type; }
	CoorType getCoorType() const { return coorType_; }
	TriCoord& operator+=(const TriCoord& pt);
	TriCoord& operator-=(const TriCoord& pt);
	TriCoord& operator/=(double scalar);
	TriCoord operator-(const TriCoord& pt)const;
	TriCoord operator+(const TriCoord& pt)const;
	bool operator<(const TriCoord& other) const;
	size_t getIndex() const { return index_; }
private:
	double x_, y_, z_;
	CoorType coorType_;
	size_t  index_;
};

struct CoorHash {
	std::size_t operator()(const TriCoord& pt) const;
};