#pragma once

#include "preprocessor.hpp"
#include "Vector3D.hpp"
#include "utils.hpp"

class Vector4D {
public:
	double x;
	double y;
	double z;
	double w;
public:
	_CExp Vector4D() : x(0), y(0), z(0), w(0) {}
	_Ex _CExp Vector4D(const Vector3D& v, double W = 0) : x(v.x), y(v.y), z(v.z), w(W) {}
	_Ex _CExp Vector4D(double X, double Y, double Z, double W = 0) : x(X), y(Y), z(Z), w(W) {}

	_NoD double length() const _NoE {
		return std::sqrt(x * x + y * y + z * z + w * w);
	}
	_NoD _CExp double lengthSqr() const _NoE {
		return x * x + y * y + z * z + w * w;
	}

	_NoD _CExp Vector4D operator+(const Vector4D& v) const _NoE {
		return Vector4D(x + v.x, y + v.y, z + v.z, w + v.w);
	}
	_NoD _CExp Vector4D operator-(const Vector4D& v) const _NoE {
		return Vector4D(x - v.x, y - v.y, z - v.z, w - v.w);
	}
	_NoD _CExp Vector4D operator*(double s) const _NoE {
		return Vector4D(x * s, y * s, z * s, w * s);
	}
	_NoD _CExp Vector4D operator/(double s) const {
		if (s == 0) {
			throw "Division by zero";
		}
		return Vector4D(x / s, y / s, z / s, w / s);
	}

	_NoD _CExp bool operator==(const Vector4D& v) const _NoE {
		return
			equalsEpsilon(v.x, x) &&
			equalsEpsilon(v.y, y) &&
			equalsEpsilon(v.z, z) &&
			equalsEpsilon(v.w, w);
	}

	_NoD _CExp Vector4D normalized() const {
		return (*this) / length();
	}
	_CExp void normalize() {
		(*this) = (*this) / length();
	}

	_NoD _CExp double dot(const Vector4D& v) const _NoE {
		return x * v.x + y * v.y + z * v.z + w * v.w;
	}
};
