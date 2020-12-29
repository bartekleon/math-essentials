#pragma once

#include "preprocessor.hpp"
#include <array>

class Matrix3x3 {
public:
	std::array<double, 9> m;
public:
	_CExp Matrix3x3() : m({
		1, 0, 0,
		0, 1, 0,
		0, 0, 1
	}) {}
	_Ex _CExp Matrix3x3(double m00, double m01, double m02, double m10, double m11, double m12, double m20, double m21, double m22) : m({
		m00, m01, m02,
		m10, m11, m12,
		m20, m21, m22
	}) { }
	_Ex _CExp Matrix3x3(const std::array<double, 9>& aflValues) : m(aflValues) {}
	_Ex _CExp Matrix3x3(std::array<double, 9>&& aflValues) : m(aflValues) {}

	_CExp void identity() _NoE {
		m[1] = m[2] = m[3] = m[5] = m[6] = m[7] = 0;
		m[0] = m[4] = m[8] = 1;
	}

	_NoD _CExp bool isIdentity() const _NoE {
		if (
			m[0] == 1 && m[1] == 0 && m[2] == 0 &&
			m[3] == 0 && m[4] == 1 && m[5] == 0 &&
			m[6] == 0 && m[7] == 0 && m[8] == 1
		) {
			return true;
		}

		return false;
	}

	_CExp void transpose() _NoE {
		(*this) = (*this).transposed();
	}
	_NoD _CExp Matrix3x3 transposed() const _NoE {
		return Matrix3x3(
			m[0], m[3], m[6],
			m[1], m[4], m[7],
			m[2], m[5], m[8]
		);
	}

	// Simple matrix operations
	_NoD _CExp Matrix3x3 operator+(const Matrix3x3& t) const _NoE {
		return Matrix3x3(
			m[0] + t.m[0], m[1] + t.m[1], m[2] + t.m[2],
			m[3] + t.m[3], m[4] + t.m[4], m[5] + t.m[5],
			m[6] + t.m[6], m[7] + t.m[7], m[8] + t.m[8]
		);
	}
	_NoD _CExp Matrix3x3 operator-(const Matrix3x3& t) const _NoE {
		return Matrix3x3(
			m[0] - t.m[0], m[1] - t.m[1], m[2] - t.m[2],
			m[3] - t.m[3], m[4] - t.m[4], m[5] - t.m[5],
			m[6] - t.m[6], m[7] - t.m[7], m[8] - t.m[8]
		);
	}
	_NoD _CExp Matrix3x3 operator*(double f) const _NoE {
		return Matrix3x3(
			m[0] * f, m[1] * f, m[2] * f,
			m[3] * f, m[4] * f, m[5] * f,
			m[6] * f, m[7] * f, m[8] * f
		);
	}
	_NoD _CExp Matrix3x3 operator/(double f) const {
		return Matrix3x3(
			m[0] / f, m[1] / f, m[2] / f,
			m[3] / f, m[4] / f, m[5] / f,
			m[6] / f, m[7] / f, m[8] / f
		);
	}

	_NoD _CExp Matrix3x3 operator*(const Matrix3x3& t) const _NoE {
		Matrix3x3 r{};

		// [a b c][A B C]   [aA+bD+cG aB+bE+cH ...
		// [d e f][D E F] = [dA+eD+fG ...
		// [g h i][G H I]

		for (int i = 0; i < 9; i += 3) {
			for (int j = 0; j < 3; j++) {
				r.m[i + j] = m[j] * t.m[i] + m[3 + j] * t.m[i + 1] + m[6 + j] * t.m[i + 2];
			}
		}

		return r;
	}

	_NoD _CExp bool operator==(const Matrix3x3& t) const _NoE {
		for (int i = 0; i < 9; i++) {
			if (!equalsEpsilon(m[i], t.m[i])) {
				return false;
			}
		}
		return true;
	}
};
