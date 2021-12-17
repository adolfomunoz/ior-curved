#pragma once

#include <../atmospheres/atmosphere.h>
#include <mj2/tracer/primitives/sphere.h>
#include <optional>
#include <list>
#include <patterns/patterns.h>

class SolverBase : public pattern::SelfRegisteringReflectableBase {
public:
    static const char* type_name() { return "solver"; }
    virtual std::list<Eigen::Vector3f> trajectory(
        const tracer::Ray& ray, const AtmosphereBase& atmosphere, const tracer::Sphere& surface) const = 0;
    virtual std::optional<tracer::Hit> trace(const tracer::Ray& ray, const AtmosphereBase& atmosphere,    const tracer::Sphere& surface) const = 0;
};

class Solver : public pattern::Pimpl<SolverBase> {
public:
    using pattern::Pimpl<SolverBase>::Pimpl;
    using pattern::Pimpl<SolverBase>::operator=;
    
    std::list<Eigen::Vector3f> trajectory(
        const tracer::Ray& ray, const AtmosphereBase& atmosphere, const tracer::Sphere& surface) const override {
            return impl()->trajectory(ray,atmosphere,surface);
        }
    virtual std::optional<tracer::Hit> trace(const tracer::Ray& ray, const AtmosphereBase& atmosphere,    const tracer::Sphere& surface) const override {
        return impl()->trace(ray,atmosphere,surface);
    }
        
};


