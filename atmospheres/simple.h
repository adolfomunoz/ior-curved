#pragma once

#include "atmosphere.h"

class Simple : public pattern::SelfRegisteringReflectable<Simple,AtmosphereBase> {
    float planet_radius;

public:
    Simple(const float planet_radius = 6370949) : planet_radius(planet_radius) {}

    float ior (float x, float y, float z) const override { // Index of Refraction
        float h = std::sqrt(x*x + y*y + z*z) - planet_radius;
        return 1.0f + 1.e-6*325*std::exp(-0.00012180 * h);
    };

    std::array<float, 3> dior (float x, float y, float z) const override { // IOR gradient
        float h = std::sqrt(x*x + y*y + z*z) - planet_radius;
        float k = 1.e-6*325*(-0.00012180)*std::exp(-0.00012180 * h)/std::sqrt(x*x+y*y+z*z);
        return std::array<float,3>{k*x,k*y,k*z};
    };
    
    static const char* type_name() { return "simple"; }

};
