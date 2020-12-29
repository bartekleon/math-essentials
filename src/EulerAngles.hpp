#pragma once

#include "preprocessor.hpp"
#include "Vector3D.hpp"
#include "coordinates.hpp"
#include "utils.hpp"

class EulerAngles {
public:
	double p;
	double y;
	double r;
public:
	_CExp EulerAngles() : p(0), y(0), r(0) {}
	_Ex _CExp EulerAngles(double pitch, double yaw, double roll) : p(pitch), y(yaw), r(roll) {}

	_NoD _CExp Vector3D toVector3D() const _NoE {
		double cp = std::cos(p);

		return Vector3D(std::cos(y) * cp, std::sin(p), std::sin(y) * cp);
	}
	_NoD _CExp Quaternion toQuaternion() const _NoE {
		return eulerAnglesToQuaternion(y, p, r);
	}

	_CExp void normalize() _NoE {
		if (p > M_PI_2) {
			p = M_PI_2;
		} else if (p < -M_PI_2) {
			p = -M_PI_2;
		}
		while (y < -M_PI) {
			y += M_PI * 2;
		}
		while (y > M_PI) {
			y -= M_PI * 2;
		}
	}
};
