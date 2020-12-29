#pragma once

#include "preprocessor.hpp"
#include "Vector3D.hpp"
#include "Quaternion.hpp"
#include "EulerAngles.hpp"

// Rotation of a vector around an axis by theta.
_NoD _CExp Vector3D rotateVectorAroundAxisAngle(const Vector3D& n, double theta, const Vector3D& v) _NoE {
    double a = std::cos(theta);
    return v * a + (n * v.dot(n) * (1 - a)) + (n.cross(v) * std::sin(theta));
}


_NoD _CExp Vector3D sphericalToCartesian(double inclination, double azimuth, double radius = 1) _NoE {
    double s = radius * std::sin(inclination);

    return Vector3D(s * std::cos(azimuth), s * std::sin(azimuth), radius * std::cos(inclination));
}


_NoD _CExp Vector3D cartesianToSpherical(double x, double y, double z) _NoE {
    double radius = std::sqrt(x * x + y * y + z * z);
    if (radius == 0) {
        return Vector3D(0, 0, 0);
    }
    double inclination = std::acos(z / radius);
    double azimuth = std::atan2(y, x);

    return Vector3D(inclination, azimuth, radius);
}
_NoD _CExp Vector3D cartesianToSpherical(const Vector3D& v) _NoE {
    return cartesianToSpherical(v.x, v.y, v.z);
}


_NoD _CExp EulerAngles quaternionToEulerAngles(double w, double x, double y, double z) _NoE {
    EulerAngles angles;
    double sinrCosp = 2 * (w * x + y * z);
    double cosrCosp = 1 - 2 * (x * x + y * y);

    // roll (x-axis rotation)
    angles.r = std::atan2(sinrCosp, cosrCosp);

    // pitch (y-axis rotation)
    double sinp = 2 * (w * y - z * x);
    if (std::abs(sinp) >= 1) {
        angles.p = std::copysign(sinp, M_PI_2);
    } else {
        angles.p = std::asin(sinp);
    }

    // yaw (z-axis rotation)
    double sinyCosp = 2 * (w * z + x * y);
    double cosyCosp = 1 - 2 * (y * y + z * z);
    angles.y = std::atan2(sinyCosp, cosyCosp);

    return angles;
}


_NoD _CExp Quaternion eulerAnglesToQuaternion(double yaw, double pitch, double roll) _NoE {
    yaw /= 2;
    pitch /= 2;
    roll /= 2;
    double cy = std::cos(yaw);
    double sy = std::sin(yaw);
    double cp = std::cos(pitch);
    double sp = std::sin(pitch);
    double cr = std::cos(roll);
    double sr = std::sin(roll);

    return Quaternion(
        cr * cp * cy + sr * sp * sy,
        sr * cp * cy - cr * sp * sy,
        cr * sp * cy + sr * cp * sy,
        cr * cp * sy - sr * sp * cy
    );
}
