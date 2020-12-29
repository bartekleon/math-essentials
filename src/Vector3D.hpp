#pragma once

#include "preprocessor.hpp"
#include "Vector2D.hpp"
#include "Matrix3x3.hpp"
#include "utils.hpp"

class Vector3D {
public:
	double x;
	double y;
	double z;
public:
	_CExp Vector3D() : x(0), y(0), z(0) {}
	_Ex _CExp Vector3D(const Vector2D& v, double Z = 0) : x(v.x), y(v.y), z(Z) {}
	_Ex _CExp Vector3D(double X, double Y, double Z) : x(X), y(Y), z(Z) {}

	_NoD double length() const _NoE {
		return std::sqrt(x * x + y * y + z * z);
	}
	_NoD _CExp double lengthSqr() const _NoE {
		return x * x + y * y + z * z;
	}

	_NoD _CExp Vector3D operator-() const _NoE {
		return Vector3D(-x, -y, -z);
	}
	_NoD _CExp Vector3D operator+(const Vector3D& v) const _NoE {
		return Vector3D(x + v.x, y + v.y, z + v.z);
	}
	_NoD _CExp Vector3D operator-(const Vector3D& v) const _NoE {
		return Vector3D(x - v.x, y - v.y, z - v.z);
	}
	_NoD _CExp Vector3D operator*(double s) const _NoE {
		return Vector3D(x * s, y * s, z * s);
	}
	_NoD _CExp Vector3D operator/(double s) const {
		if (s == 0) {
			throw "Division by zero";
		}
		return Vector3D(x / s, y / s, z / s);
	}
	_CExp Vector3D& operator+=(const Vector3D& v) _NoE {
		x += v.x;
		y += v.y;
		z += v.z;
		return *this;
	}
	_CExp Vector3D& operator-=(const Vector3D& v) _NoE {
		x -= v.x;
		y -= v.y;
		z -= v.z;
		return *this;
	}
	_CExp Vector3D& operator*=(double s) _NoE {
		x *= s;
		y *= s;
		z *= s;
		return *this;
	}
	_CExp Vector3D& operator/=(double s) {
		if (s == 0) {
			throw "Division by zero";
		}
		x /= s;
		y /= s;
		z /= s;
		return *this;
	}
	_CExp Vector3D& operator*(const Vector3D& v) const _NoE {
		return Vector3D(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
	}
	_NoD _CExp bool operator==(const Vector3D& v) const _NoE {
		return
			equalsEpsilon(x, v.x) &&
			equalsEpsilon(y, v.y) &&
			equalsEpsilon(z, v.z);
	}

	_NoD _CExp Vector3D normalized() const {
		return (*this) / length();
	}
	_CExp void normalize() {
		(*this) = (*this) / length();
	}

	_NoD _CExp double dot(const Vector3D& v) const _NoE {
		return x * v.x + y * v.y + z * v.z;
	}
	_NoD _CExp Vector3D cross(const Vector3D& v) const _NoE {
		return Vector3D(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
	}
	_NoD _CExp Matrix3x3 outer(const Vector3D& v) const _NoE {
		return Matrix3x3(
			x * v.x, x * v.y, x * v.z,
			y * v.x, y * v.y, y * v.z,
			z * v.x, z * v.y, z * v.z
		);
	}
};
