#pragma once

#include "gtest/gtest.h"

#include "../src/Vector3D.hpp"

TEST(Vector3D, ShouldCalculateLength) {
	const Vector3D v{ 2.4231, 34.123, 6.132 };

	EXPECT_DOUBLE_EQ(v.length(), std::sqrt(2.4231 * 2.4231 + 34.123 * 34.123 + 6.132 * 6.132));
	EXPECT_DOUBLE_EQ(v.lengthSqr(), 2.4231 * 2.4231 + 34.123 * 34.123 + 6.132 * 6.132);
}

TEST(Vector3D, ShouldCheckIfVectorsAreEqual) {
	const Vector3D v1{ 2.4231, 34.123, 1.2345 };
	const Vector3D v2{ -2.4231, 4.123, 324.123 };
	const Vector3D v3{ 2.42310000, 34.1230000, 1.23450000 };

	EXPECT_TRUE(v1 == v3);
	EXPECT_FALSE(v1 == v2);
}

TEST(Vector3D, ShouldPerformBasicOperations) {
	const Vector3D v1{ 2.4231, 34.123, 123.543 };
	const Vector3D v2{ -2.4231, 4.123, -0.4216 };
	const Vector3D v3{ 2.42310000, 34.1230000, 123.543 };

	EXPECT_TRUE(v1 + v2 == Vector3D(0, 38.246, 123.1214));
	EXPECT_TRUE(v1 - v3 == Vector3D{});
	EXPECT_TRUE(v2 * 3 == Vector3D(-7.2693, 12.369, -1.2648));
	EXPECT_TRUE(v1 / 2 == Vector3D(1.21155, 17.0615, 61.7715));
	EXPECT_ANY_THROW(v1 / 0);
}

TEST(Vector3D, ShouldNormalize) {
	const Vector3D v1{ 2.4231, 34.123, 43.13354 };
	const Vector3D v2{ 0, 0, 0 };

	EXPECT_DOUBLE_EQ(v1.normalized().length(), 1);
	EXPECT_ANY_THROW(v2.normalized());
}

TEST(Vector3D, ShouldPerformDotProduct) {
	const Vector3D v{ 2.4231, 34.123, 1.23 };

	EXPECT_DOUBLE_EQ(v.dot(v), 1171.76344261);
}

TEST(Vector3D, ShouldPerformCrossProduct) {
	const Vector3D v1{ 2.4231, 34.123, 1.23 };
	const Vector3D v2{ -0.342, 1.342, 1.23 };

	EXPECT_TRUE(v1.cross(v1) == Vector3D{});
	EXPECT_TRUE(v1.cross(v2) == Vector3D(40.32063, -3.401073, 14.9218662));
}

TEST(Vector3D, ShouldPerformOuterProduct) {
	const Vector3D v{ 2.4231, 34.123, 1.23 };
	const Matrix3x3 expected{
		5.87141361, 82.6834413, 2.980413,
		82.6834413, 1164.379129, 41.97129,
		2.980413, 41.97129, 1.5129
	};
	EXPECT_TRUE(v.outer(v) == expected);
}
