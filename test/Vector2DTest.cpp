#pragma once

#include "gtest/gtest.h"

#include "../src/Vector2D.hpp"

TEST(Vector2D, ShouldCalculateLength) {
	const Vector2D v{ 2.4231, 34.123 };

	EXPECT_DOUBLE_EQ(v.length(), std::sqrt(2.4231 * 2.4231 + 34.123 * 34.123));
	EXPECT_DOUBLE_EQ(v.lengthSqr(), 2.4231 * 2.4231 + 34.123 * 34.123);
}

TEST(Vector2D, ShouldCheckIfVectorsAreEqual) {
	const Vector2D v1{ 2.4231, 34.123 };
	const Vector2D v2{ -2.4231, 4.123 };
	const Vector2D v3{ 2.42310000, 34.1230000 };

	EXPECT_TRUE(v1 == v3);
	EXPECT_FALSE(v1 == v2);
}

TEST(Vector2D, ShouldPerformBasicOperations) {
	const Vector2D v1{ 2.4231, 34.123 };
	const Vector2D v2{ -2.4231, 4.123 };
	const Vector2D v3{ 2.42310000, 34.1230000 };

	EXPECT_TRUE(v1 + v2 == Vector2D(0, 34.123 + 4.123));
	EXPECT_TRUE(v1 - v3 == Vector2D{});
	EXPECT_TRUE(v2 * 3 == Vector2D(-7.2693, 12.369));
	EXPECT_TRUE(v1 / 2 == Vector2D(1.21155, 17.0615));
	EXPECT_ANY_THROW(v1 / 0);
}

TEST(Vector2D, ShouldNormalize) {
	const Vector2D v1{ 2.4231, 34.123 };
	const Vector2D v2{ 0, 0 };

	EXPECT_DOUBLE_EQ(v1.normalized().length(), 1);
	EXPECT_ANY_THROW(v2.normalized());
}

TEST(Vector2D, ShouldPerformDotProduct) {
	const Vector2D v{ 2.4231, 34.123 };

	EXPECT_DOUBLE_EQ(v.dot(v), 1170.25054261);
}
