#include <svg-cpp-plot/svg-cpp-plot.h>
#include <ivp/ivp.h>
#include <mj2/tracer/primitives/sphere.h>
#include "../eq/fermat.h"
#include <cmath>
#include <string>

#include "../atmospheres/simple.h"
#include "../atmospheres/strange.h"

int main(int argc, char** argv) {
    /** Initialization */
    /** Avaiable parameters:
     * -zenith <angle> // Zenith angle
     * -omega_ch <omega_ch> // Cherenkov angle
     */

    float width = 200; // Plot width
    float height = 200; // Plot height
    const double radius = 6370949; // Earth radius (m)
    const float detector_height = 0; // Meters adobe radius
    const float method_steps = 1000.0f; // Steps for the numeric method

    auto ior = [=] (float x, float y, float z) { // Index of Refraction
        auto atmosphere = Simple(radius);
        return atmosphere.ior(x,y,z);
    };
    auto dior = [=] (float x, float y, float z) { // IOR gradient
        auto atmosphere = Simple(radius);
        return atmosphere.dior(x,y,z);
    };

    // Numeric Method for fermat solving
    auto function = fermat(ior,dior);
    IVP::RungeKutta2 method(method_steps);

    // Zenith
    float angle = 88*M_PI/180.0; // Default
    for (int i = 0; i<(argc-1); ++i) { // Custom angle as argument
        if (std::string(argv[i]) == "-zenith") angle = atof(argv[++i])*M_PI/180.0;
    }

    // Scene
    tracer::Sphere earth(Eigen::Vector3f(0,0,0),radius);
    tracer::Sphere surface(Eigen::Vector3f(0,0,0),radius+detector_height);

    // Get atmosphere border (origin of the ray) from the Earth with the zenith angle
    float distance = (-23.3f + 19.45f/std::cos(0.99273f*angle)) * 1000; // x1000 km -> m
    float x = distance * std::sin(angle);
    float y = 0;
    float z = ( distance * std::cos(angle) ) + radius + detector_height;
    float atmosphere_height = std::sqrt(x*x + y*y + z*z) - radius;
    Eigen::Vector3f origin(x,0,z);

    // Cherenkov cone. For -1, 0 and 1 angles
    float omega_ch = std::acos(1/ior(x,0,z))*M_PI/180.0; // Default
    for (int i = 0; i<(argc-1); ++i) { // Custom angle as argument
        if (std::string(argv[i]) == "-omega_ch") omega_ch = atof(argv[++i])*M_PI/180.0;
    }

    // Print initial information
    printf("Zenith: %f degrees\n", angle*(180/M_PI));
    printf("Cherenkov angle: %f degrees\n", omega_ch*(180/M_PI));
    printf("Origin of rays: %f,%f\n", x, z);
    printf("Distance Origin-North: %f m\n", distance);
    printf("Atmosphere height: %f m\n", atmosphere_height);






    /** Ray tracing (& path plotting) **/
    svg_cpp_plot::SVGPlot plt;
    std::list<float> hits_x, hits_y, nohits_x, nohits_y, hits_nonlinear_x, hits_nonlinear_y;
    for (float a = (angle-omega_ch); a<=(angle+1.5*omega_ch); a+=omega_ch) {
        std::cout << "Angle used: " << a << "\n";

        // Trace ray towards the Earth
        tracer::Ray ray(origin, Eigen::Vector3f(-std::sin(a), 0, -std::cos(a)));
        if (auto hit = surface.trace(ray)) {
            hits_x.push_back((*hit).point()[0]);
            hits_y.push_back((*hit).point()[2]);
        } else {
            auto point = ray.at(radius);
            nohits_x.push_back(point[0]);
            nohits_y.push_back(point[2]);
        }

        // Trace curved ray towards de Earth
        Eigen::Array<float, 6, 1> ini;
        std::list<float> path_x, path_y;
        ini(Eigen::seq(Eigen::fix<0>, Eigen::fix<2>)) = ray.origin();
        ini(Eigen::seq(Eigen::fix<3>, Eigen::fix<5>)) = ray.direction();
        for (auto s : method.steps(function, 0.0f, ini, 4.0f * distance)) {
            ray.set_range_max(s.step());
            if (auto hit = surface.trace(ray)) { // Hit
                hits_nonlinear_x.push_back((*hit).point()[0]);
                hits_nonlinear_y.push_back((*hit).point()[2]);
                path_x.push_back(hits_nonlinear_x.back());
                path_y.push_back(hits_nonlinear_y.back());
                break;
            } else { // No hit
                path_x.push_back(s.y()[0]);
                path_y.push_back(s.y()[2]);
                if (std::abs(s.y()[1]) > 1.e-2) std::cerr << "Warning : displacement in y :" << s.y()[1] << std::endl;
                ray = tracer::Ray(s.y()(Eigen::seq(Eigen::fix<0>, Eigen::fix<2>)),
                                  s.y()(Eigen::seq(Eigen::fix<3>, Eigen::fix<5>)));
            }
        }

        // Plot path of curved ray from origin towards Earth
        plt.plot(path_x,path_y).color("r").linewidth(0.25);
        //plt.scatter(path_x,path_y).c("r").s(0.2);
    }

    /** Plotting **/
    plt.scatter(hits_x,hits_y).c("k"); // Plot point where ray intersect Earth
    plt.scatter(hits_nonlinear_x,hits_nonlinear_y).c("r"); // Plot point where curved ray intersect Earth
    plt.scatter({x}, {z}).c("g"); // Plot origin of the rays

    // Plot Earth
    std::list<float> earth_x, earth_y;
    for (float a = -std::asin(x/radius); a<=std::asin(x/radius); a+=M_PI/180) {
        earth_x.push_back(radius*std::sin(a));
        earth_y.push_back(radius*std::cos(a));
    }
    plt.plot(earth_x,earth_y).color("b").linewidth(0.1);

    // Plot Surface
    std::list<float> surface_x, surface_y;
    float surface_height = 3000; // m
    for (float a = -std::asin(x/(radius+surface_height)); a<=std::asin(x/(radius+surface_height)); a+=M_PI/180) {
        surface_x.push_back((radius+surface_height)*std::sin(a));
        surface_y.push_back((radius+surface_height)*std::cos(a));
    }
    plt.plot(surface_x,surface_y).color("b").linewidth(0.1);

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

    // Vertical line at x=0
    plt.plot({0,0},{-10e8,10e8}).color("k").linewidth(0.5).alpha(0.25);

    // Create figure and save it.
    printf("Creating figure...");
    plt.axis(limits).linewidth(0).xticks({}).yticks({}).figsize({width,height}).savefig("atmosphere.svg");

}
