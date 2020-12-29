#pragma once

#include "preprocessor.hpp"
#include "utils.hpp"

#include <cmath>

class Vector2D {
public:
	double x;
	double y;
public:
	_CExp Vector2D() : x(0), y(0) {}
	_Ex _CExp Vector2D(double X, double Y) : x(X), y(Y) {}

	_NoD double length() const _NoE {
		return std::sqrt(x * x + y * y);
	}
	_NoD _CExp double lengthSqr() const _NoE {
		return x * x + y * y;
	}

	_NoD _CExp Vector2D operator-() const _NoE {
		return Vector2D(-x, -y);
	}
	_NoD _CExp Vector2D operator+(const Vector2D& v) const _NoE {
		return Vector2D(x + v.x, y + v.y);
	}
	_NoD _CExp Vector2D operator-(const Vector2D& v) const _NoE {
		return Vector2D(x - v.x, y - v.y);
	}
	_NoD _CExp Vector2D operator*(double s) const _NoE {
		return Vector2D(x * s, y * s);
	}
	_NoD _CExp Vector2D operator/(double s) const {
		if (s == 0) {
			throw "Division by zero";
		}
		return Vector2D(x / s, y / s);
	}
	_CExp Vector2D& operator+=(const Vector2D& v) _NoE {
		x += v.x;
		y += v.y;
		return *this;
	}
	_CExp Vector2D& operator-=(const Vector2D& v) _NoE {
		x -= v.x;
		y -= v.y;
		return *this;
	}
	_CExp Vector2D& operator*=(double s) _NoE {
		x *= s;
		y *= s;
		return *this;
	}
	_CExp Vector2D& operator/=(double s) {
		if (s == 0) {
			throw "Division by zero";
		}
		x /= s;
		y /= s;
		return *this;
	}

	_NoD _CExp bool operator==(const Vector2D& v) const _NoE {
		return
			equalsEpsilon(x, v.x) &&
			equalsEpsilon(y, v.y);
	}

	_NoD _CExp Vector2D normalized() const {
		return (*this) / length();
	}
	_CExp void normalize() {
		(*this) = (*this) / length();
	}

	_NoD _CExp double dot(const Vector2D& v) const _NoE {
		return x * v.x + y * v.y;
	}
};
