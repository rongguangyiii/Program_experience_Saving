#pragma once
#include <vector>
#include <iomanip>

struct Coord
{
	double x_, y_, z_;
	Coord(double x, double y, double z = 0.0) : x_(x), y_(y), z_(z){}
	Coord() : x_(0), y_(0), z_(0) {}

public:
	const double& x() const { return x_; }
	const double& y() const { return y_; }
	const double& z() const { return z_; }
	Coord& operator/=(double scalar);
	Coord& operator+=(const Coord& other);
	Coord& operator-=(const Coord& other);
	bool operator==(const Coord& other) const;
	bool operator<(const Coord& other) const;
	void operator=(const Coord& other);

};
