#pragma once

#include "gtest/gtest.h"

#include "../src/Vector4D.hpp"

TEST(Vector4D, ShouldCalculateLength) {
	const Vector4D v{ 2.4231, 34.123, 6.132, 1.23 };

	double result = 2.4231 * 2.4231 + 34.123 * 34.123 + 6.132 * 6.132 + 1.23 * 1.23;

	EXPECT_DOUBLE_EQ(v.length(), std::sqrt(result));
	EXPECT_DOUBLE_EQ(v.lengthSqr(), result);
}

TEST(Vector4D, ShouldCheckIfVectorsAreEqual) {
	const Vector4D v1{ 2.4231, 34.123, 1.2345, 0.1234 };
	const Vector4D v2{ -2.4231, 4.123, 324.123, 0.43123 };
	const Vector4D v3{ 2.42310000, 34.1230000, 1.23450000, 0.1234 };

	EXPECT_TRUE(v1 == v3);
	EXPECT_FALSE(v1 == v2);
}

TEST(Vector4D, ShouldPerformBasicOperations) {
	const Vector4D v1{ 2.4231, 34.123, 123.543, -0.512 };
	const Vector4D v2{ -2.4231, 4.123, -0.4216, 1.232 };
	const Vector4D v3{ 2.42310000, 34.1230000, 123.543, -0.512 };

	EXPECT_TRUE(v1 + v2 == Vector4D(0, 38.246, 123.1214, 0.72));
	EXPECT_TRUE(v1 - v3 == Vector4D{});
	EXPECT_TRUE(v2 * 3 == Vector4D(-7.2693, 12.369, -1.2648, 3.696));
	EXPECT_TRUE(v1 / 2 == Vector4D(1.21155, 17.0615, 61.7715, -0.256));
	EXPECT_ANY_THROW(v1 / 0);
}

TEST(Vector4D, ShouldNormalize) {
	const Vector4D v1{ 2.4231, 34.123, 43.13354, 0.342 };
	const Vector4D v2{ 0, 0, 0, 0 };

	EXPECT_DOUBLE_EQ(v1.normalized().length(), 1);
	EXPECT_ANY_THROW(v2.normalized());
}
