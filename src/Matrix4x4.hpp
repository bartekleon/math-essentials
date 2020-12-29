#pragma once

#include "preprocessor.hpp"
#include "Vector3D.hpp"
#include "Vector4D.hpp"

#include <array>

class Matrix4x4 {
public:
	std::array<double, 16> m;
public:
	_CExp Matrix4x4() : m({
		1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1
	}) {}
	_Ex _CExp Matrix4x4(double m00, double m01, double m02, double m03, double m10, double m11, double m12, double m13, double m20, double m21, double m22, double m23, double m30, double m31, double m32, double m33) : m({
		m00, m01, m02, m03,
		m10, m11, m12, m13,
		m20, m21, m22, m23,
		m30, m31, m32, m33
	}) { }
	_Ex _CExp Matrix4x4(const std::array<double, 16>& aflValues) : m(aflValues) {}
	_Ex _CExp Matrix4x4(std::array<double, 16>&& aflValues) : m(aflValues) {}
	_Ex _CExp Matrix4x4(const Vector3D& vecForward, const Vector3D& vecUp, const Vector3D& vecRight, const Vector3D& vecPosition = Vector3D(0, 0, 0)) : m({
		vecForward.x, vecForward.y, vecForward.z, 0,
		vecUp.x, vecUp.y, vecUp.z, 0,
		vecRight.x, vecRight.y, vecRight.z, 0,
		vecPosition.x, vecPosition.y, vecPosition.z, 1
	}) { }

	_CExp void identity() _NoE {
		m.fill(0);
		m[0] = m[5] = m[10] = m[15] = 1;
	}

	_NoD _CExp bool isIdentity() const _NoE {
		if (
			m[0] == 1 && m[1] == 0 && m[2] == 0 && m[3] == 0 &&
			m[4] == 0 && m[5] == 1 && m[6] == 0 && m[7] == 0 &&
			m[8] == 0 && m[9] == 0 && m[10] == 1 && m[11] == 0 &&
			m[12] == 0 && m[13] == 0 && m[14] == 0 && m[15] == 1
		) {
			return true;
		}

		return false;
	}

	_CExp void transpose() _NoE {
		(*this) = (*this).transposed();
	}
	_NoD _CExp Matrix4x4 transposed() const _NoE {
		return Matrix4x4(
			m[0], m[4], m[8], m[12],
			m[1], m[5], m[9], m[13],
			m[2], m[6], m[10], m[14],
			m[3], m[7], m[11], m[15]
		);
	}

	// Simple matrix operations
	_NoD _CExp Matrix4x4 operator+(const Matrix4x4& t) const _NoE {
		return Matrix4x4(
			m[0] + t.m[0], m[1] + t.m[1], m[2] + t.m[2], m[3] + t.m[3],
			m[4] + t.m[4], m[5] + t.m[5], m[6] + t.m[6], m[7] + t.m[7],
			m[8] + t.m[8], m[9] + t.m[9], m[10] + t.m[10], m[11] + t.m[11],
			m[12] + t.m[12], m[13] + t.m[13], m[14] + t.m[14], m[15] + t.m[15]
		);
	}
	_NoD _CExp Matrix4x4 operator-(const Matrix4x4& t) const _NoE {
		return Matrix4x4(
			m[0] - t.m[0], m[1] - t.m[1], m[2] - t.m[2], m[3] - t.m[3],
			m[4] - t.m[4], m[5] - t.m[5], m[6] - t.m[6], m[7] - t.m[7],
			m[8] - t.m[8], m[9] - t.m[9], m[10] - t.m[10], m[11] - t.m[11],
			m[12] - t.m[12], m[13] - t.m[13], m[14] - t.m[14], m[15] - t.m[15]
		);
	}
	_NoD _CExp Matrix4x4 operator*(double f) const _NoE {
		return Matrix4x4(
			m[0] * f, m[1] * f, m[2] * f, m[3] * f,
			m[4] * f, m[5] * f, m[6] * f, m[7] * f,
			m[8] * f, m[9] * f, m[10] * f, m[11] * f,
			m[12] * f, m[13] * f, m[14] * f, m[15] * f
		);
	}
	_NoD _CExp Matrix4x4 operator/(double f) const {
		return Matrix4x4(
			m[0] / f, m[1] / f, m[2] / f, m[3] / f,
			m[4] / f, m[5] / f, m[6] / f, m[7] / f,
			m[8] / f, m[9] / f, m[10] / f, m[11] / f,
			m[12] / f, m[13] / f, m[14] / f, m[15] / f
		);
	}

	// Set a transformation
	_CExp void setTranslation(const Vector3D& vecPos) _NoE {
		m[12] = vecPos.x;
		m[13] = vecPos.y;
		m[14] = vecPos.z;
	}
	_CExp void setRotation(double radians, const Vector3D& vecAxis) _NoE {
		const Vector3D v = vecAxis.normalized();

		// c = cos(angle), s = sin(angle), t = (1-c)
		// [ xxt+c   xyt-zs  xzt+ys ]
		// [ yxt+zs  yyt+c   yzt-xs ]
		// [ zxt-ys  zyt+xs  zzt+c  ]

		double x = v.x;
		double y = v.y;
		double z = v.z;

		double c = std::cos(radians);
		double s = std::sin(radians);
		double t = 1 - c;

		m[0] = x * x * t + c;
		m[4] = x * y * t - z * s;
		m[8] = x * z * t + y * s;

		m[1] = y * x * t + z * s;
		m[5] = y * y * t + c;
		m[9] = y * z * t - x * s;

		m[2] = z * x * t - y * s;
		m[6] = z * y * t + x * s;
		m[10] = z * z * t + c;
	}
	_CExp void setScale(const Vector3D& vecScale) _NoE {
		m[0] = vecScale.x;
		m[5] = vecScale.y;
		m[10] = vecScale.z;
	}
	// Reflection around a plane with this normal which passes through the center of the local space.
	_CExp void setReflection(const Vector3D& vecPlaneNormal) _NoE {
		const Vector3D vecPlane = vecPlaneNormal.normalized();

		m[0] = 1 - 2 * vecPlane.x * vecPlane.x;
		m[5] = 1 - 2 * vecPlane.y * vecPlane.y;
		m[10] = 1 - 2 * vecPlane.z * vecPlane.z;
		m[4] = m[1] = -2 * vecPlane.x * vecPlane.y;
		m[8] = m[2] = -2 * vecPlane.x * vecPlane.z;
		m[6] = m[9] = -2 * vecPlane.y * vecPlane.z;
	}

	// Just like gluPerspectives
	_NoD _CExp static Matrix4x4 projectPerspective(double flFOV, double flAspectRatio, double flNear, double flFar) {
		double flRight = flNear * std::tan(flFOV / 2);
		double flLeft = -flRight;

		double flBottom = flLeft / flAspectRatio;
		double flTop = flRight / flAspectRatio;

		return projectFrustum(flLeft, flRight, flBottom, flTop, flNear, flFar);
	}
	// Just like glFrustum
	_NoD _CExp static Matrix4x4 projectFrustum(double flLeft, double flRight, double flBottom, double flTop, double flNear, double flFar) {
		Matrix4x4 m{};

		double flXD = flRight - flLeft;
		double flYD = flTop - flBottom;
		double flZD = flFar - flNear;

		m.m[0] = (2 * flNear) / flXD;
		m.m[5] = (2 * flNear) / flYD;

		m.m[8] = (flRight + flLeft) / flXD;
		m.m[9] = (flTop + flBottom) / flYD;
		m.m[10] = -(flFar + flNear) / flZD;
		m.m[11] = -1;

		m.m[14] = -(2 * flFar * flNear) / flZD;

		m.m[15] = 0;

		return m;
	}
	// Just like glOrtho
	_NoD _CExp static Matrix4x4 projectOrthographic(double flLeft, double flRight, double flBottom, double flTop, double flNear, double flFar) {
		Matrix4x4 m{};

		double flXD = flRight - flLeft;
		double flYD = flTop - flBottom;
		double flZD = flFar - flNear;

		m.m[0] = 2.0f / flXD;
		m.m[5] = 2.0f / flYD;
		m.m[10] = -2.0f / flZD;

		m.m[12] = -(flRight + flLeft) / flXD;
		m.m[13] = -(flTop + flBottom) / flYD;
		m.m[14] = -(flFar + flNear) / flZD;

		return m;
	}
	// Like gluLookAt but a direction parameter instead of target
	_NoD _CExp static Matrix4x4 constructCameraView(const Vector3D& vecPosition, const Vector3D& vecDirection, const Vector3D& vecUp) {
		Matrix4x4 m{};

		const Vector3D vecDir = vecDirection.normalized();

		const Vector3D vecCamSide = vecDir.cross(vecUp).normalized();
		const Vector3D vecCamUp = vecCamSide.cross(vecDir);

		m.setForwardVector(Vector3D(vecCamSide.x, vecCamUp.x, -vecDir.x));
		m.setUpVector(Vector3D(vecCamSide.y, vecCamUp.y, -vecDir.y));
		m.setRightVector(Vector3D(vecCamSide.z, vecCamUp.z, -vecDir.z));

		m.addTranslation(-vecPosition);

		return m;
	}

	// Add a translation
	_CExp Matrix4x4 operator+=(const Vector3D& v) _NoE {
		m[12] += v.x;
		m[13] += v.y;
		m[14] += v.z;

		return *this;
	}
	_CExp Matrix4x4 operator+(const Vector3D& v) const _NoE {
		Matrix4x4 r = *this;

		r.m[12] += v.x;
		r.m[13] += v.y;
		r.m[14] += v.z;

		return r;
	}
	_CExp Matrix4x4 operator-=(const Vector3D& v) _NoE {
		m[12] -= v.x;
		m[13] -= v.y;
		m[14] -= v.z;

		return *this;
	}
	_CExp Matrix4x4 operator-(const Vector3D& v) const _NoE {
		Matrix4x4 r = *this;

		r.m[12] -= v.x;
		r.m[13] -= v.y;
		r.m[14] -= v.z;

		return r;
	}

	// Add a transformation
	_NoD _CExp Matrix4x4 operator*(const Matrix4x4& t) const _NoE {
		Matrix4x4 r{};

		// [a b c d][A B C D]   [aA+bE+cI+dM
		// [e f g h][E F G H] = [eA+fE+gI+hM ...
		// [i j k l][I J K L]
		// [m n o p][M N O P]

		for (int i = 0; i < 16; i += 4)	{
			for (int j = 0; j < 4; j++) {
				r.m[i + j] = m[j] * t.m[i] + m[4 + j] * t.m[i + 1] + m[8 + j] * t.m[i + 2] + m[12 + j] * t.m[i + 3];
			}
		}

		return r;
	}
	_NoD _CExp Matrix4x4 operator*=(const Matrix4x4& t) _NoE {
		*this = (*this) * t;

		return *this;
	}

	_NoD _CExp bool operator==(const Matrix4x4& t) const _NoE {
		for (int i = 0; i < 16; i++) {
			if (std::fabs(m[i] - t.m[i]) > std::numeric_limits<double>::epsilon()) {
				return false;
			}
		}

		return true;
	}
	_NoD _CExp bool equals(const Matrix4x4& t, double epsilon = std::numeric_limits<double>::epsilon()) const _NoE {
		for (int i = 0; i < 16; i++) {
			if (std::fabs(m[i] - t.m[i]) > std::numeric_limits<double>::epsilon()) {
				return false;
			}
		}

		return true;
	}

	// Add a transformation
	_NoD _CExp Matrix4x4 addTranslation(const Vector3D& v) _NoE {
		Matrix4x4 r{};
		r.setTranslation(v);
		(*this) *= r;

		return *this;
	}
	_NoD _CExp Matrix4x4 addScale(const Vector3D& vecScale) _NoE {
		Matrix4x4 r{};
		r.setScale(vecScale);
		(*this) *= r;

		return *this;
	}
	_NoD _CExp Matrix4x4 addReflection(const Vector3D& vecPlaneNormal) _NoE {
		Matrix4x4 r{};
		r.setReflection(vecPlaneNormal);
		(*this) *= r;

		return *this;
	}

	_NoD _CExp Vector3D getTranslation() const _NoE {
		return Vector3D(m[12], m[13], m[14]);
	}
	_NoD _CExp Vector3D getScale() const _NoE {
		return Vector3D(getForwardVector().length(), getUpVector().length(), getRightVector().length());
	}

	// Transform a vector
	_NoD _CExp Vector3D operator*(const Vector3D& v) const _NoE {
		// [a b c x][X] 
		// [d e f y][Y] = [aX+bY+cZ+x dX+eY+fZ+y gX+hY+iZ+z]
		// [g h i z][Z]
		//          [1]
		Vector3D vecResult;
		vecResult.x = m[0] * v.x + m[4] * v.y + m[8] * v.z + m[12];
		vecResult.y = m[1] * v.x + m[5] * v.y + m[9] * v.z + m[13];
		vecResult.z = m[2] * v.x + m[6] * v.y + m[10] * v.z + m[14];
		return vecResult;
	}
	// Same as homogenous vector with w=0 transform, no translation.
	// You want to use this for directional vectors such as normals and velocities because translations will change their length.
	// It's not immune to scaling though! A matrix with scaling will output a vector of different length than the input.
	_NoD _CExp Vector3D transformVector(const Vector3D& v) const _NoE {
		// [a b c][X] 
		// [d e f][Y] = [aX+bY+cZ dX+eY+fZ gX+hY+iZ]
		// [g h i][Z]
		Vector3D vecResult;
		vecResult.x = m[0] * v.x + m[4] * v.y + m[8] * v.z;
		vecResult.y = m[1] * v.x + m[5] * v.y + m[9] * v.z;
		vecResult.z = m[2] * v.x + m[6] * v.y + m[10] * v.z;
		return vecResult;
	}

	_NoD _CExp Vector4D operator*(const Vector4D& v) const _NoE {
		// [a b c x][X] 
		// [d e f y][Y] = [aX+bY+cZ+xW dX+eY+fZ+yW gX+hY+iZ+zW jX+kY+lZ+mW]
		// [g h i z][Z]
		// [j k l m][W]
		Vector4D vecResult;
		vecResult.x = m[0] * v.x + m[4] * v.y + m[8] * v.z + m[12] * v.w;
		vecResult.y = m[1] * v.x + m[5] * v.y + m[9] * v.z + m[13] * v.w;
		vecResult.z = m[2] * v.x + m[6] * v.y + m[10] * v.z + m[14] * v.w;
		vecResult.w = m[3] * v.x + m[7] * v.y + m[11] * v.z + m[15] * v.w;
		return vecResult;
	}

	// Try not to use these in case the underlying format changes.
	_NoD _CExp Vector4D getRow(int i) const {
		i *= 4;
		return Vector4D(m[i], m[i + 1], m[i + 2], m[i + 3]);
	}
	_NoD _CExp Vector4D getColumn(int i) const {
		return Vector4D(m[i], m[4 + i], m[8 + i], m[12 + i]);
	}
	_CExp void setColumn(int i, const Vector4D& vecColumn) {
		m[i] = vecColumn.x;
		m[4 + i] = vecColumn.y;
		m[8 + i] = vecColumn.z;
		m[12 + i] = vecColumn.w;
	}
	_CExp void setColumn(int i, const Vector3D& vecColumn) {
		m[i] = vecColumn.x;
		m[4 + i] = vecColumn.y;
		m[8 + i] = vecColumn.z;
	}

	_CExp void setForwardVector(const Vector3D& vecForward) _NoE {
		m[0] = vecForward.x;
		m[1] = vecForward.y;
		m[2] = vecForward.z;
	}

	_CExp void setUpVector(const Vector3D& vecUp) _NoE {
		m[4] = vecUp.x;
		m[5] = vecUp.y;
		m[6] = vecUp.z;
	}
	_CExp void setRightVector(const Vector3D& vecRight) _NoE {
		m[8] = vecRight.x;
		m[9] = vecRight.y;
		m[10] = vecRight.z;
	}
	_NoD _CExp Vector3D getForwardVector() const _NoE  {
		return Vector3D(m[0], m[1], m[2]);
	}
	_NoD _CExp Vector3D getUpVector() const _NoE {
		return Vector3D(m[4], m[5], m[6]);
	}
	_NoD _CExp Vector3D getRightVector() const _NoE {
		return Vector3D(m[8], m[9], m[10]);
	}

	// Not a true inversion, only works if the matrix is a translation/rotation matrix.
	_CExp void invertRT() _NoE {
		Matrix4x4 t{};

		for (int h = 0; h < 3; h++) {
			for (int v = 0; v < 3; v++) {
				t.m[4 * h + v] = m[4 * v + h];
			}
		}

		Vector3D vecTranslation = -getTranslation();

		m = t.m;

		setTranslation(t * vecTranslation);
	}
	_NoD _CExp Matrix4x4 invertedRT() const _NoE {
		Matrix4x4 r{};

		for (int h = 0; h < 3; h++)
			for (int v = 0; v < 3; v++)
				r.m[4 * h + v] = m[4 * v + h];

		r.setTranslation(r * (-getTranslation()));

		return r;
	}

	_NoD _CExp double trace() const _NoE {
		return m[0] + m[5] + m[9];
	}
};
