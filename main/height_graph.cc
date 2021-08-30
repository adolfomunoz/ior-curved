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
    const float method_steps = 100.0f; // Steps for the numeric method

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






    /** Ray tracing */
    std::list<float> heights, x_points, diors;
    std::cout << "angle: " << angle << "\n";
    tracer::Ray ray(origin, Eigen::Vector3f(-std::sin(angle), 0, -std::cos(angle)));
    Eigen::Array<float, 6, 1> ini;
    ini(Eigen::seq(Eigen::fix<0>, Eigen::fix<2>)) = ray.origin();
    ini(Eigen::seq(Eigen::fix<3>, Eigen::fix<5>)) = ray.direction();

    // Trace curved ray towards de Earth
    for (auto s : method.steps(function, 0.0f, ini, 4.0f * distance)) {//1.0f*float(atmosphere_height))) {
        ray.set_range_max(s.step());
        if (auto hit = surface.trace(ray)) {
            float height = std::sqrt((*hit).point()[0]*(*hit).point()[0] + (*hit).point()[1]*(*hit).point()[1] + (*hit).point()[2]*(*hit).point()[2]) - radius;
            heights.push_back(height);
            x_points.push_back((*hit).point()[0]);
            diors.push_back(dior((*hit).point()[0], (*hit).point()[1], (*hit).point()[2])[0]);
            break;
        } else {
            float height = std::sqrt(s.y()[0]*s.y()[0] + s.y()[1]*s.y()[1] + s.y()[2]*s.y()[2]) - radius;
            heights.push_back(height);
            x_points.push_back(s.y()[0]);
            diors.push_back(dior(s.y()[0], s.y()[1], s.y()[2])[0]);
            if (std::abs(s.y()[1]) > 1.e-2) std::cerr << "Warning : displacement in y :" << s.y()[1] << std::endl;
            ray = tracer::Ray(s.y()(Eigen::seq(Eigen::fix<0>, Eigen::fix<2>)),
                              s.y()(Eigen::seq(Eigen::fix<3>, Eigen::fix<5>)));
        }
    }

    /** Plotting Height graph **/
    svg_cpp_plot::SVGPlot plt_height;
    plt_height.scatter(x_points,heights);

    // Markers for different sections on atmosphere definition
    /*plt_height.plot({0,333000},{600,600}).color("k").linewidth(0.5).alpha(0.25);
    plt_height.plot({0,333000},{425,425}).color("k").linewidth(0.5).alpha(0.25);
    plt_height.plot({0,333000},{350,350}).color("k").linewidth(0.5).alpha(0.25);
    plt_height.plot({0,333000},{275,275}).color("k").linewidth(0.5).alpha(0.25);*/

    // Create figure and save it.
    printf("Creating figure (Height)...");
    plt_height.xlabel("Distance to pole (rho or x-axis)");
    plt_height.ylabel("Height from sea level");
    plt_height.title("Height of the ray coming to the Earth");
    plt_height.axis({0,40000,0,1000}).figsize({width,height}).savefig("height_graph.svg");


    /** Plotting dIOR graph **/
    svg_cpp_plot::SVGPlot plt_dior;
    plt_dior.scatter(heights,diors);

    // Create figure and save it.
    printf("Creating figure (dIOR)...");
    plt_height.xlabel("Height from the sea level");
    plt_height.ylabel("dIOR");
    plt_height.title("dIOR for a ray coming to the Earth");
    plt_dior.axis({0,25000, -1e-9,1e-9}).figsize({width,height}).savefig("dior_graph.svg");


}
