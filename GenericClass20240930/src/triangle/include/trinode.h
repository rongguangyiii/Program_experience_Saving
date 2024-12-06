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
#ifndef TRINODE_H_
#define TRINODE_H_
#include<algorithm>
#include<cmath>
#include<stdexcept>

enum class TriNodeType
{
	Notset,
	outter,
	inner,
};

class TriNode
{
public:
	TriNode() :x_(0), y_(0), z_(0), index_(0), type_(TriNodeType::Notset) {};
	TriNode(double x, double y, double z = 0.0, size_t id = 0, TriNodeType type = TriNodeType::Notset) :x_(x), y_(y), z_(z), index_(id), type_(type) {}
	TriNode(const TriNode& pt) :x_(pt.x_), y_(pt.y_), z_(pt.z_), type_(pt.type_), index_(pt.index_) {}
	~TriNode() {};
	void setTriNodeIdx(size_t id) { index_ = id; }
	void setTriNodeType(TriNodeType type) { type_ = type; }
	TriNodeType gettype()const { return type_; }
	const size_t getindex() const { return index_; }
public:
	const double& x() const { return x_; }
	const double& y() const { return y_; }
	const double& z() const { return z_; }
	void SetCoord(const double& x = 0, const double& y = 0, const double& z = 0);
	bool operator ==(const TriNode& pt)const;
	double norm() const;
	TriNode& operator+=(const TriNode& pt);
	TriNode& operator-=(const TriNode& pt);
	TriNode& operator/=(double scalar);
	TriNode operator-(const TriNode& pt)const;
	TriNode operator+(const TriNode& pt)const;
	bool operator<(const TriNode& other) const;
private:
	size_t index_;
	double x_, y_, z_;
	TriNodeType type_;

};

#endif // !CDTPOINT_H_
