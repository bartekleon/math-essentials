#pragma once

#include "preprocessor.hpp"

#include <limits>

_CExp long double PI_RAD =  0.0174532925199432957L;
_CExp long double PI_DEG = 57.2957795130823208L;

_NoD _CExp double deg_to_rad(double degrees) _NoE {
    return degrees * PI_RAD;
}

_NoD _CExp double rad_to_deg(double rad) _NoE {
    return rad * PI_DEG;
}

_CExp double epsilon = std::numeric_limits<double>::epsilon();

_NoD _CExp bool equalsEpsilon(double v1, double v2, double eps = epsilon) _NoE {
    double r = v1 - v2;
    if (r < 0) {
        r = -r;
    }
    if (v1 < 0) {
        v1 = -v1;
    }
    if (v2 < 0) {
        v2 = -v2;
    }
    return r <= (v1 > v2 ? v2 : v1) * eps;
}
