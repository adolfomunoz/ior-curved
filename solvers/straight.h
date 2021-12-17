#pragma once

#include "solver.h"

class Straight : public pattern::SelfRegisteringReflectable<Straight,SolverBase> {
 
public:
    Straight() {} //An empty constructor is required, otherwise it does not work
    
    std::list<Eigen::Vector3f> trajectory(
        const tracer::Ray& ray, const AtmosphereBase& atmosphere, const tracer::Sphere& surface) const override {
            std::list<Eigen::Vector3f> t;
            t.push_back(ray.origin());
            if (auto hit = surface.trace(ray)) t.push_back((*hit).point());
            return t;
        }
        
    std::optional<tracer::Hit> trace(
        const tracer::Ray& ray, const AtmosphereBase& atmosphere, const tracer::Sphere& surface) const override {
            return surface.trace(ray);
        }
        
    static const char* type_name() { return "straight"; }

};
