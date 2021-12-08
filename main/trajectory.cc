#include <svg-cpp-plot/svg-cpp-plot.h>
#include <ivp/ivp.h>
#include <mj2/tracer/primitives/sphere.h>
#include "../eq/fermat.h"
#include <cmath>
#include <string>

#include "../atmospheres/simple.h"
#include "../atmospheres/strange.h"
#include "../atmospheres/gutierrez2006.h"
#include "../solvers/straight.h"
#include "../solvers/ivp.h"

/**
 * Image that shows (side face) the Earth and the path of straight and curved rays towards the surface.
 */
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

    Atmosphere atmosphere = Simple();
    pattern::load_commandline(atmosphere,argc,argv);
    auto solvers = pattern::make_from_commandline<std::list<Solver>>(argc,argv);

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
    float omega_ch = std::acos(1/atmosphere.ior(x,0,z))*M_PI/180.0; // Default
    for (int i = 0; i<(argc-1); ++i) { // Custom angle as argument
        if (std::string(argv[i]) == "-omega_ch") omega_ch = atof(argv[++i])*M_PI/180.0;
    }

    // Print initial information
    printf("Zenith: %f degrees\n", angle*(180/M_PI));
    printf("Cherenkov angle: %f degrees\n", omega_ch*(180/M_PI));
    printf("Origin of rays: %f,%f\n", x, z);
    printf("Distance Origin-North: %f m\n", distance);
    printf("Atmosphere height: %f m\n", atmosphere_height);


    std::vector<std::string> colors{"k","r","g","b"};



    /** Ray tracing (& path plotting) **/
    svg_cpp_plot::SVGPlot plt; int c = 0;
    float min_x = x; float min_y = z;
    for (auto solver : solvers) {
        const std::string& color = colors[c++];
        for (float a = (angle-omega_ch); a<=(angle+1.5*omega_ch); a+=omega_ch) {
            // Trace ray towards the Earth
            tracer::Ray ray(origin, Eigen::Vector3f(-std::sin(a), 0, -std::cos(a)));
            auto trajectory = solver.trajectory(ray,atmosphere,surface);
            auto hit = solver.trace(ray,atmosphere,surface);
            std::list<float> xs, ys;
            for (auto point : trajectory) { xs.push_back(point[0]); ys.push_back(point[2]); }
            plt.plot(xs,ys).color(color).linewidth(0.25);
            if (hit) {  
                plt.scatter({float((*hit).point()[0])},{float((*hit).point()[2])}).c(color);
                if ((*hit).point()[0]<min_x) min_x = (*hit).point()[0];
                if ((*hit).point()[2]<min_y) min_y = (*hit).point()[2];
            }                
        }    
    }
    
    //Plot origin
    plt.scatter({origin[0]},{origin[2]}).c("k").marker("s");

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
    std::array<float,4> limits{min_x,x,min_y,z};
    if ((limits[1]-limits[0])*width>(limits[3]-limits[2])*height)
        limits[2] = limits[3] - (limits[1]-limits[0])*height/width;
    else
        limits[0] = limits[1] - (limits[3]-limits[2])*width/height;

    float d = (limits[3]-limits[2])/32.0;
    limits[0] -= d; limits[1] += d; limits[2] -=d; limits[3] += d;

    // Vertical line at x=0
    plt.plot({0,0},{-10e8,10e8}).color("k").linewidth(0.5).alpha(0.25);

    // Create figure and save it.
    printf("Creating figure...");
    plt.axis(limits).linewidth(0).xticks({}).yticks({}).figsize({width,height}).savefig("atmosphere.svg");

}
