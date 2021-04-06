#include <svg-cpp-plot/svg-cpp-plot.h>
#include <ivp/ivp.h>
#include <mj2/tracer/primitives/sphere.h>
#include "../eq/fermat.h"
#include <cmath>

int main() {
    const double radius = 6370949;
    const double angle = 88*M_PI/180.0;
    const double atmosphere_height = 325000;
    float width = 200;
    float height = 100;
    tracer::Sphere earth(Eigen::Vector3f(0,0,0),radius);
    auto ior = [=] (float x, float y, float z) {
        float h = std::sqrt(x*x + y*y + z*z) - radius;
        return 1.0f + 1.e-6*325*std::exp(-0.00012180 * h);
    };
    
    auto dior = [=] (float x, float y, float z) {
        float h = std::sqrt(x*x + y*y + z*z) - radius;
        float k = 1.e-6*325*(-0.00012180)*std::exp(-0.00012180 * h)/std::sqrt(x*x+y*y+z*z);
        return std::array<float,3>{k*x,k*y,k*z};
    };
    
    //To solve where the origin of the ray is (x,z coordinates) we need to solve an order 2 equation and keep the positive root
    double tana = std::tan(angle);
    double tana2 = tana*tana;
    
    double a = tana2 + 1;
    double b = -2.0f*radius*tana2;
    double c = radius*radius*tana2 - (radius+height)*(radius+atmosphere_height);

    float z = (-b + std::sqrt(b*b - 4.0f*a*c))/(2*a); 
    float x = (z-radius)*tana;
    
    Eigen::Vector3f origin(x,0,z);
    
    auto function = fermat(ior,dior);
    IVP::Adaptive<IVP::Dopri> method(100,1.e-3);
//    IVP::Euler method(100);
   
    svg_cpp_plot::SVGPlot plt;
    std::list<float> hits_x, hits_y, nohits_x, nohits_y, hits_nonlinear_x, hits_nonlinear_y;
    for (float a = (angle-M_PI/180.0); a<=(angle+1.5*M_PI/180.0); a+=M_PI/180.0) {
        tracer::Ray ray(origin,Eigen::Vector3f(-std::sin(a),0,-std::cos(a)));
        if (auto hit = earth.trace(ray)) {
            hits_x.push_back((*hit).point()[0]);
            hits_y.push_back((*hit).point()[2]);
        } else {
            auto point = ray.at(radius);
            nohits_x.push_back(point[0]);
            nohits_y.push_back(point[2]);
        }
        std::list<float> path_x, path_y;
        Eigen::Array<float,6,1> ini; 
        ini(Eigen::seq(Eigen::fix<0>,Eigen::fix<2>)) = ray.origin();
        ini(Eigen::seq(Eigen::fix<3>,Eigen::fix<5>)) = ray.direction();
        for (auto s : method.steps(function,0.0f,ini,10.0f*float(atmosphere_height))) {
            ray.set_range_max(s.step());
            if (auto hit = earth.trace(ray)) {
                hits_nonlinear_x.push_back((*hit).point()[0]);
                hits_nonlinear_y.push_back((*hit).point()[2]);
                path_x.push_back(hits_nonlinear_x.back());
                path_y.push_back(hits_nonlinear_y.back());
                break;
            } else {
                path_x.push_back(s.y()[0]);
                path_y.push_back(s.y()[2]);
                ray = tracer::Ray(s.y()(Eigen::seq(Eigen::fix<0>,Eigen::fix<2>)),s.y()(Eigen::seq(Eigen::fix<3>,Eigen::fix<5>)));
            }
        }
        plt.plot(path_x,path_y).color("r").linewidth(0.25);
    }
    
    std::list<float> earth_x, earth_y;
    for (float a = -std::asin(x/radius); a<=std::asin(x/radius); a+=M_PI/180) {
        earth_x.push_back(radius*std::sin(a));
        earth_y.push_back(radius*std::cos(a));
    }
    
    float min_x = *(std::min_element(hits_x.begin(),hits_x.end()));
    float min_y = *(std::min_element(hits_y.begin(),hits_y.end()));
    std::array<float,4> limits{min_x,x,min_y,z};
    if ((limits[1]-limits[0])*width>(limits[3]-limits[2])*height)
        limits[2] = limits[3] - (limits[1]-limits[0])*height/width;
    else
        limits[0] = limits[1] - (limits[3]-limits[2])*width/height;
        
    float d = (limits[3]-limits[2])/32.0;
    limits[0] -= d; limits[1] += d; limits[2] -=d; limits[3] += d;
    
    plt.plot(earth_x,earth_y).color("b");
    
    std::list<float>::const_iterator i,j;
    for (i = hits_x.begin(), j = hits_y.begin(); (i != hits_x.end()) && (j != hits_y.end()); ++i, ++j)
        plt.plot({x,*i},{z,*j}).color("k").linewidth(0.25).format("--");
    for (i = nohits_x.begin(), j = nohits_y.begin(); (i != nohits_x.end()) && (j != nohits_y.end()); ++i, ++j)
        plt.plot({x,*i},{z,*j}).color("k").linewidth(0.25).format("--");
        
    plt.scatter(hits_x,hits_y).c("k");
    plt.scatter(hits_nonlinear_x,hits_nonlinear_y).c("r");
    plt.scatter({x}, {z}).c("g");;
    
    
    
    plt.axis(limits).linewidth(0).xticks({}).yticks({}).figsize({width,height}).savefig("atmosphere.svg");

}