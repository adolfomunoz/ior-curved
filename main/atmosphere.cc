#include <svg-cpp-plot/svg-cpp-plot.h>
#include <ivp/ivp.h>
#include <mj2/tracer/primitives/sphere.h>
#include "../eq/fermat.h"
#include <cmath>
#include <string>
// Tested values: -zenith 88 -height 325000 -omega_ch 1
int main(int argc, char** argv) {
    /** Parameters */
    // Plot
    float width = 200;
    float height =200;

    // Rays
    double angle = 88*M_PI/180.0; // Zenith angle to radians
    for (int i = 0; i<(argc-1); ++i) { // Custom angle as argument
        if (std::string(argv[i]) == "-zenith") angle = atof(argv[++i])*M_PI/180.0;
    }
    printf("Zenith: %f degrees\n", angle*(180/M_PI));

    // Scene
    const double radius = 6370949; // Earth radius
    tracer::Sphere earth(Eigen::Vector3f(0,0,0),radius);

    auto ior = [=] (float x, float y, float z) { // Index of Refraction
        float h = std::sqrt(x*x + y*y + z*z) - radius; // h = zv
        return 1.0f + 1.e-6*325*std::exp(-0.00012180 * h);
    };
    auto dior = [=] (float x, float y, float z) { // IOR gradient
        float h = std::sqrt(x*x + y*y + z*z) - radius;
        float k = 1.e-6*325*(-0.00012180)*std::exp(-0.00012180 * h)/std::sqrt(x*x+y*y+z*z);
        return std::array<float,3>{k*x,k*y,k*z};
    };

    /** Get atmosphere border (origin of the ray) from the Earth with the zenith angle */
    //To solve where the origin of the ray is (x,z coordinates) we need to solve an order 2 equation and keep the positive root
    float distance = (-10 + 12 / std::cos(angle)) * 1000; // x1000 km -> m
    float x = distance * std::sin(angle);
    float z = ( distance * std::cos(angle) ) + radius;
    float atmosphere_height = std::sqrt(x*x + z*z) - radius;

    printf("Origin of rays: %f,%f\n", x, z);
    Eigen::Vector3f origin(x,0,z);
    printf("Distance Origin-North: %f m\n", distance);
    printf("Atmosphere height: %f m\n", atmosphere_height);

    auto function = fermat(ior,dior);
    //IVP::Adaptive<IVP::Dopri> method(100,1.e-3); // (Pasos inicialmente, tolerancia)
    IVP::RungeKutta2 method(2000.0f);

    /** Cherenkov cone */
    // For -1, 0 and 1 angles
    double omega_ch = 88*M_PI/180.0; //std::acos(1/ior(x,0,z))*M_PI/180.0; //arccos(1/ior)
    for (int i = 0; i<(argc-1); ++i) { // Custom angle as argument
        if (std::string(argv[i]) == "-omega_ch") omega_ch = atof(argv[++i])*M_PI/180.0;
    }
    printf("omega_ch: %f degrees\n", omega_ch*(180/M_PI));

    /** Plotting */
    svg_cpp_plot::SVGPlot plt;
    std::list<float> hits_x, hits_y, nohits_x, nohits_y, hits_nonlinear_x, hits_nonlinear_y;
    for (float a = (angle-omega_ch); a<=(angle+1.5*omega_ch); a+=omega_ch) {
        //for (float a = angle; a <= angle*1.5; a += angle) {
        tracer::Ray ray(origin, Eigen::Vector3f(-std::sin(a), 0, -std::cos(a)));
        // Trace ray towards the Earth
        if (auto hit = earth.trace(ray)) {
            hits_x.push_back((*hit).point()[0]);
            hits_y.push_back((*hit).point()[2]);
        } else {
            auto point = ray.at(radius);
            nohits_x.push_back(point[0]);
            nohits_y.push_back(point[2]);
        }
        std::list<float> path_x, path_y;
        Eigen::Array<float, 6, 1> ini;
        ini(Eigen::seq(Eigen::fix<0>, Eigen::fix<2>)) = ray.origin();
        ini(Eigen::seq(Eigen::fix<3>, Eigen::fix<5>)) = ray.direction();
        // Trace curved ray towards de Earth
        for (auto s : method.steps(function, 0.0f, ini, 4.0f * distance)) {//1.0f*float(atmosphere_height))) {
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
                if (std::abs(s.y()[1]) > 1.e-2) std::cerr << "Warning : displacement in y :" << s.y()[1] << std::endl;
                ray = tracer::Ray(s.y()(Eigen::seq(Eigen::fix<0>, Eigen::fix<2>)),
                                  s.y()(Eigen::seq(Eigen::fix<3>, Eigen::fix<5>)));
            }
        }
        // Plot curved ray from origin towards Earth
        plt.plot(path_x,path_y).color("r").linewidth(0.25);
        //plt.scatter(path_x,path_y).c("r").s(0.2);
    }

    // Plot Earth
    std::list<float> earth_x, earth_y;
    for (float a = -std::asin(x/radius); a<=std::asin(x/radius); a+=M_PI/180) {
        earth_x.push_back(radius*std::sin(a));
        earth_y.push_back(radius*std::cos(a));
    }
    plt.plot(earth_x,earth_y).color("b").linewidth(0.1);

    // Get the limits of the plot
    float min_x = *(std::min_element(hits_x.begin(),hits_x.end()));
    float min_y = *(std::min_element(hits_y.begin(),hits_y.end()));
    std::array<float,4> limits{min_x,x,min_y,z};
    if ((limits[1]-limits[0])*width>(limits[3]-limits[2])*height)
        limits[2] = limits[3] - (limits[1]-limits[0])*height/width;
    else
        limits[0] = limits[1] - (limits[3]-limits[2])*width/height;

    float d = (limits[3]-limits[2])/32.0;
    limits[0] -= d; limits[1] += d; limits[2] -=d; limits[3] += d;


    // Plot rays from origin towards Earth as discontinued line
    std::list<float>::const_iterator i,j;
    for (i = hits_x.begin(), j = hits_y.begin(); (i != hits_x.end()) && (j != hits_y.end()); ++i, ++j)
        plt.plot({x,*i},{z,*j}).color("k").linewidth(0.25).format("--");
    for (i = nohits_x.begin(), j = nohits_y.begin(); (i != nohits_x.end()) && (j != nohits_y.end()); ++i, ++j)
        plt.plot({x,*i},{z,*j}).color("k").linewidth(0.25).format("--");

    plt.scatter(hits_x,hits_y).c("k"); // Plot point where ray intersect Earth
    plt.scatter(hits_nonlinear_x,hits_nonlinear_y).c("r"); // Plot point where curved ray intersect Earth
    plt.scatter({x}, {z}).c("g"); // Plot origin of the rays


    // Create figure and save it.
    printf("Creating figure...");
    plt.axis(limits).linewidth(0).xticks({}).yticks({}).figsize({width,height}).savefig("atmosphere.svg");

}
