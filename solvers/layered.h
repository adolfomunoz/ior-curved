#pragma once

#include "solver.h"

class LayeredSolver : public pattern::SelfRegisteringReflectable<LayeredSolver,SolverBase> {
    float step;
    unsigned long max_steps;
public:
    LayeredSolver(float step = 0.1, unsigned long max_steps = 1000000) :
        step(step), max_steps(max_steps) {}
    auto reflect() { return std::tie(step,max_steps); }
    auto reflect_names() const { return std::tuple("step","max-steps"); }
    
    std::list<Eigen::Vector3f> trajectory(
        const tracer::Ray& r, const AtmosphereBase& atmosphere, const tracer::Sphere& surface) const override {
            std::list<Eigen::Vector3f> t;
            tracer::Ray ray = r;
            
            for (unsigned long s = 0; s<max_steps; ++s) {
                float distance = ray.origin().norm();
                t.push_back(ray.origin());
                float ior = 1.0f;
                Eigen::Vector3f layerhit;
                Eigen::Vector3f layernormal;
                tracer::Sphere nextlayer(Eigen::Vector3f(0,0,0),distance+step);
                tracer::Sphere lastlayer(Eigen::Vector3f(0,0,0),distance+step);
                if (auto hit = tracer::Sphere(Eigen::Vector3f(0,0,0),distance-step).trace(ray)) {
                    ray.set_range_max((*hit).distance());
                    layerhit = (*hit).point();
                    layernormal=(*hit).normal();
                    Eigen::Vector3f dhit = layerhit/layerhit.norm();
                    Eigen::Vector3f queryior1 = layerhit + 0.5*step*dhit;
                    Eigen::Vector3f queryior2 = layerhit - 0.5*step*dhit;
                    ior = atmosphere.ior(queryior1[0],queryior1[1],queryior1[2])/
                            atmosphere.ior(queryior2[0],queryior2[1],queryior2[2]);
                }
                if (auto hit = tracer::Sphere(Eigen::Vector3f(0,0,0),distance+step).trace(ray)) {
                    ray.set_range_max((*hit).distance());
                    layerhit = (*hit).point();
                    layernormal=(*hit).normal(); 
                    Eigen::Vector3f dhit = layerhit/layerhit.norm();
                    Eigen::Vector3f queryior1 = layerhit - 0.5*step*dhit;
                    Eigen::Vector3f queryior2 = layerhit + 0.5*step*dhit;
                    ior = atmosphere.ior(queryior1[0],queryior1[1],queryior1[2])/
                            atmosphere.ior(queryior2[0],queryior2[1],queryior2[2]);
                }
                if (auto hit = surface.trace(ray)) {
                    t.push_back((*hit).point()); break;
                } else {
                    float disc = 1.0f - ior*ior*layernormal.cross(ray.direction()).squaredNorm();
                    if (disc >= 0.0) 
                        ray = tracer::Ray(layerhit,ior*layernormal.cross(-layernormal.cross(ray.direction()))-layernormal*std::sqrt(disc));
                    else //Total internal reflection
                        ray = tracer::Ray(layerhit,ray.direction()-2.0f*ray.direction().dot(layernormal)*layernormal);
                }
               
            }
            return t;
        }
        
    std::optional<tracer::Hit> trace(
        const tracer::Ray& r, const AtmosphereBase& atmosphere, const tracer::Sphere& surface) const override {
            tracer::Ray ray = r;
            
            for (unsigned long s = 0; s<max_steps; ++s) {
                float distance = ray.origin().norm();          
                float ior = 1.0f;
                Eigen::Vector3f layerhit;
                Eigen::Vector3f layernormal;
                tracer::Sphere nextlayer(Eigen::Vector3f(0,0,0),distance+step);
                tracer::Sphere lastlayer(Eigen::Vector3f(0,0,0),distance+step);
                if (auto hit = tracer::Sphere(Eigen::Vector3f(0,0,0),distance-step).trace(ray)) {
                    ray.set_range_max((*hit).distance());
                    layerhit = (*hit).point();
                    layernormal=(*hit).normal();
                    Eigen::Vector3f dhit = layerhit/layerhit.norm();
                    Eigen::Vector3f queryior1 = layerhit + 0.5*step*dhit;
                    Eigen::Vector3f queryior2 = layerhit - 0.5*step*dhit;
                    ior = atmosphere.ior(queryior1[0],queryior1[1],queryior1[2])/
                            atmosphere.ior(queryior2[0],queryior2[1],queryior2[2]);
                }
                if (auto hit = tracer::Sphere(Eigen::Vector3f(0,0,0),distance+step).trace(ray)) {
                    ray.set_range_max((*hit).distance());
                    layerhit = (*hit).point();
                    layernormal=(*hit).normal();                    
                    Eigen::Vector3f dhit = layerhit/layerhit.norm();
                    Eigen::Vector3f queryior1 = layerhit - 0.5*step*dhit;
                    Eigen::Vector3f queryior2 = layerhit + 0.5*step*dhit;
                    ior = atmosphere.ior(queryior1[0],queryior1[1],queryior1[2])/
                            atmosphere.ior(queryior2[0],queryior2[1],queryior2[2]);
                }
                if (auto hit = surface.trace(ray)) {
                    return hit;
                } else {
                    float disc = 1.0f - ior*ior*layernormal.cross(ray.direction()).squaredNorm();
                    if (disc >= 0.0) 
                        ray = tracer::Ray(layerhit,ior*layernormal.cross(-layernormal.cross(ray.direction()))-layernormal*std::sqrt(disc));
                    else //Total internal reflection
                        ray = tracer::Ray(layerhit,ray.direction()-2.0f*ray.direction().dot(layernormal)*layernormal);
                }
               
            }
            return {};
        }
        
    static const char* type_name() { return "layered"; }
};
