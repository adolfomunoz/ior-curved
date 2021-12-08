#pragma once

#include <ivp/ivp.h>
#include "solver.h"

namespace {
    template<typename Method>
    struct method_traits { };
    
    template<>
    struct method_traits<IVP::Euler> { static const char* name() { return "euler"; } }; 
    
    template<>
    struct method_traits<IVP::RungeKutta2> { static const char* name() { return "rungekutta2"; } }; 
    
    template<>
    struct method_traits<IVP::RungeKutta4> { static const char* name() { return "rungekutta4"; } }; 
}

template<typename Method>
class IVPSolver : public pattern::SelfRegisteringReflectable<IVPSolver<Method>,SolverBase> {
    float step;
    unsigned long max_steps;
public:
    IVPSolver(float step = 0.1, unsigned long max_steps = 1000000) :
        step(step), max_steps(max_steps) {}
    auto reflect() { return std::tie(step,max_steps); }
    auto reflect_names() const { return std::tuple("step","max-steps"); }
    
    std::list<Eigen::Vector3f> trajectory(
        const tracer::Ray& r, const AtmosphereBase& atmosphere, const tracer::Sphere& surface) const override {
            std::list<Eigen::Vector3f> t;
            auto function = fermat(
                [&atmosphere] (float x, float y, float z) { return atmosphere.ior(x,y,z);},
                [&atmosphere] (float x, float y, float z) { return atmosphere.dior(x,y,z);}
            );
            tracer::Ray ray = r;
            Method method(step);
            Eigen::Array<float, 6, 1> ini;
            ini(Eigen::seq(Eigen::fix<0>, Eigen::fix<2>)) = ray.origin();
            ini(Eigen::seq(Eigen::fix<3>, Eigen::fix<5>)) = ray.direction();
            for (auto s : method.steps(function, 0.0f, ini, step*float(max_steps))) {
                t.push_back(ray.origin());
                ray.set_range_max(s.step());
                if (auto hit = surface.trace(ray)) { t.push_back((*hit).point()); break;}
                else {
                    ray = tracer::Ray(
                            s.y()(Eigen::seq(Eigen::fix<0>, Eigen::fix<2>)),
                            s.y()(Eigen::seq(Eigen::fix<3>, Eigen::fix<5>)));
                }
            }
            return t;
        }
        
    std::optional<tracer::Hit> trace(
        const tracer::Ray& r, const AtmosphereBase& atmosphere, const tracer::Sphere& surface) const override {
            auto function = fermat(
                [&atmosphere] (float x, float y, float z) { return atmosphere.ior(x,y,z);},
                [&atmosphere] (float x, float y, float z) { return atmosphere.dior(x,y,z);}
            );
            tracer::Ray ray = r;
            Method method(step);
            Eigen::Array<float, 6, 1> ini;
            ini(Eigen::seq(Eigen::fix<0>, Eigen::fix<2>)) = ray.origin();
            ini(Eigen::seq(Eigen::fix<3>, Eigen::fix<5>)) = ray.direction();
            for (auto s : method.steps(function, 0.0f, ini, step*float(max_steps))) {
                ray.set_range_max(s.step());
                if (auto hit = surface.trace(ray)) { return hit; }
                else {
                    ray = tracer::Ray(
                            s.y()(Eigen::seq(Eigen::fix<0>, Eigen::fix<2>)),
                            s.y()(Eigen::seq(Eigen::fix<3>, Eigen::fix<5>)));
                }
            }
            return {};
        }
        
    static const char* type_name() { return method_traits<Method>::name(); }
};

template class IVPSolver<IVP::Euler>;
template class IVPSolver<IVP::RungeKutta2>;
template class IVPSolver<IVP::RungeKutta4>;      
