#pragma once

#include "preprocessor.hpp"
#include "Vector2D.hpp"
#include "Vector3D.hpp"
#include "utils.hpp"

class Point {
public:
	double x;
	double y;
	double z;
public:
	_CExp Point() : x(0), y(0), z(0) {}
	_Ex _CExp Point(double X, double Y, double Z = 0) : x(X), y(X), z(Z) {}
	_Ex _CExp Point(const Vector2D& v) : x(v.x), y(v.y), z(0) {}
	_Ex _CExp Point(const Vector3D& v) : x(v.x), y(v.y), z(v.z) {}

	_NoD _CExp Point operator+(const Vector3D& v) const _NoE {
		return Point(x + v.x, y + v.y, z + v.z);
	}

	_NoD _CExp Point operator-(const Vector3D& v) const _NoE {
		return Point(x - v.x, y - v.y, z - v.z);
	}

	_NoD _CExp Point operator+(const Vector2D& v) const _NoE {
		return Point(x + v.x, y + v.y, z);
	}

	_NoD _CExp Point operator-(const Vector2D& v) const _NoE {
		return Point(x - v.x, y - v.y, z);
	}

	_NoD _CExp bool operator==(const Point& p) const _NoE {
		return
			equalsEpsilon(x, p.x) &&
			equalsEpsilon(y, p.y) &&
			equalsEpsilon(z, p.z);
	}
};
