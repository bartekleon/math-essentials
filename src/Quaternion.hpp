#pragma once

#include "preprocessor.hpp"
#include "Vector3D.hpp"
#include "Vector4D.hpp"
#include "coordinates.hpp"

class Quaternion {
public:
	double w = 1;
	double x = 0;
	double y = 0;
	double z = 0;
public:
	_CExp Quaternion() {}
	_Ex _CExp Quaternion(double W, double X, double Y, double Z) : w(W), x(X), y(Y), z(Z) {}
	_Ex _CExp Quaternion(const Vector3D& n, double rad) {
		rad /= 2;
		double s = std::sin(rad);

		w = cos(rad);
		x = n.x * s;
		y = n.y * s;
		z = n.z * s;
	}

	_NoD _CExp Quaternion operator+(const Quaternion& q) const _NoE {
		return Quaternion(w + q.w, x + q.x, y + q.y, z + q.z);
	}
	_NoD _CExp Quaternion operator-(const Quaternion& q) const _NoE {
		return Quaternion(w - q.w, x - q.x, y - q.y, z - q.z);
	}
	_NoD _CExp Quaternion operator*(double s) const _NoE {
		return Quaternion(w * s, x * s, y * s, z * s);
	}
	_NoD _CExp Quaternion operator/(double s) const {
		return Quaternion(w / s, x / s, y / s, z / s);
	}
	_NoD _CExp Quaternion operator*(const Quaternion& q) const _NoE {
		const Vector3D v1{ x, y, z };
		const Vector3D v2{ q.x, q.y, q.z };
		const Vector3D v3 = v2 * w + v1 * q.w + v1 * v2;

		return Quaternion(w * q.w - v1.dot(v2), v3.x, v3.y, v3.z);
	}
	_NoD _CExp double dot(const Quaternion& q) const _NoE {
		return Vector3D(x, y, z).dot(Vector3D(q.x, q.y, q.z)) + w * q.w;
	}

	_CExp Vector3D rotate(const Vector3D& p) const {
		double s = 2 / lengthSqr();
		return Vector3D(
			(1 - s * (y * y + z * z)) * p.x + s * (x * y - z * w) * p.y + s * (x * z + y * w) * p.z,
			s * (x * y + z * w) * p.x + (1 - s * (x * x + z * z)) * p.y + s * (y * z - x * w) * p.z,
			s * (x * z - y * w) * p.x + s * (y * z + x * w) * p.y + (1 - s * (x * x + y * y)) * p.z
		);
	}

	_NoD _CExp Quaternion normalized() const {
		return (*this) / Vector3D(x, y, z).dot(Vector3D(x, y, z));
	}
	_CExp void normalize() {
		(*this) = (*this) / Vector3D(x, y, z).dot(Vector3D(x, y, z));
	}

	_NoD _CExp Quaternion inversed() const _NoE {
		Quaternion q = *this;
		q.inverse();
		return q;
	}
	_CExp void inverse() _NoE {
		x = -x;
		y = -y;
		z = -z;
	}

	_NoD _CExp double length() const _NoE {
		return std::sqrt(w * w + x * x + y * y + z * z);
	}
	_NoD _CExp double lengthSqr() const _NoE {
		return w * w + x * x + y * y + z * z;
	}

	_NoD _CExp bool isUnit(double epsi = epsilon) const _NoE {
		return std::fabs(lengthSqr() - 1) < epsi;
	}

	_NoD _CExp EulerAngles toEulerAngles() const _NoE {
		return quaternionToEulerAngles(w, x, y, z);
	}

	_NoD _CExp Vector4D toAxisAngle() const {
		double l = std::sqrt(x * x + y * y + z * z);
		return Vector4D(x / l, y / l, z / l, 2 * std::atan2(l, w));
	}
};
