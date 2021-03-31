#include <svg-cpp-plot/svg-cpp-plot.h>
#include <ivp/ivp.h>
#include <mj2/tracer/primitives/sphere.h>
#include <cmath>

int main() {
    const double radius = 6370949;
    const double angle = 88*M_PI/180.0;
    const double atmosphere_height = 325000;
    float width = 200;
    float height = 100;
    tracer::Sphere earth(Eigen::Vector3f(0,0,0),radius);
    
    //To solve where the origin of the ray is (x,z coordinates) we need to solve an order 2 equation and keep the positive root
    double tana = std::tan(angle);
    double tana2 = tana*tana;
    
    double a = tana2 + 1;
    double b = -2.0f*radius*tana2;
    double c = radius*radius*tana2 - (radius+height)*(radius+atmosphere_height);

    float z = (-b + std::sqrt(b*b - 4.0f*a*c))/(2*a); 
    float x = (z-radius)*tana;
    
    Eigen::Vector3f origin(x,0,z);
    
   
    svg_cpp_plot::SVGPlot plt;
    std::list<float> hits_x, hits_y, nohits_x, nohits_y;
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
    
    plt.plot(earth_x,earth_y);
    
    std::list<float>::const_iterator i,j;
    for (i = hits_x.begin(), j = hits_y.begin(); (i != hits_x.end()) && (j != hits_y.end()); ++i, ++j)
        plt.plot({x,*i},{z,*j}).color("k").linewidth(0.25);
    for (i = nohits_x.begin(), j = nohits_y.begin(); (i != nohits_x.end()) && (j != nohits_y.end()); ++i, ++j)
        plt.plot({x,*i},{z,*j}).color("k").linewidth(0.25);
    
    plt.scatter(hits_x,hits_y);
    plt.scatter({x}, {z});
    
    plt.axis(limits).linewidth(0).xticks({}).yticks({}).figsize({width,height}).savefig("atmosphere.svg");
//    plt.scatter({x},{-z});
//    plt.savefig("atmosphere2.svg");

}