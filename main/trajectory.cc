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

struct Config : public pattern::Reflectable<Config> {
    float plot_width = 200;
    float plot_height = 200;
    float detector_height = 0;
    float zenith = 88; //In degrees, change to radians
    float omega_ch = -1; //Negative means generate automatic value
    std::string output = "atmosphere.svg";
    
    auto reflect() { return std::tie(plot_width,plot_height,detector_height,zenith,omega_ch,output); }
    auto reflect_names() const { return std::tuple("plot-width","plot_height","detector-height","zenith","omega-ch","output"); }
    static const char* type_name() { return "config"; } 
};

/**
 * Image that shows (side face) the Earth and the path of straight and curved rays towards the surface.
 */
int main(int argc, char** argv) {
    /** Initialization */
    /** Avaiable parameters:
     * -zenith <angle> // Zenith angle
     * -omega_ch <omega_ch> // Cherenkov angle
     */
    std::cout<<argv[0]<<std::endl;
    std::cout<<"Available atmospheres : "<<Atmosphere::registered()<<std::endl;
    std::cout<<"Available solvers     : "<<Solver::registered()<<std::endl;
    


    Atmosphere atmosphere = Simple();
    pattern::load_commandline(atmosphere,argc,argv);
    auto solvers = pattern::make_from_commandline<std::list<Solver>>(argc,argv);
    Config config;
    pattern::load_commandline(config,argc,argv);

    const double radius = 6370949; // Earth radius (m)
    float detector_height = config.detector_height; // Meters adobe radius

    // Zenith
    float angle = config.zenith*M_PI/180.0; // To radians

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
    float omega_ch = ((config.omega_ch>0)?config.omega_ch:std::acos(1/atmosphere.ior(x,0,z)))*M_PI/180.0; // If it has been set, use the set value, or use a different value otherwise

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
    
     std::cerr<<"ATMOSPHERE"<<std::endl<<atmosphere.xml()<<std::endl;
    for (auto solver : solvers) {
        std::cout<<"SOLVER"<<std::endl<<solver.xml()<<std::endl;
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

    // Plot spheres
    std::list<float> earth_x, earth_y, detector_x, detector_y, atmosphere_x, atmosphere_y;
    for (float a = -std::asin(2*x/radius); a<=std::asin(2*x/radius); a+=(std::asin(2*x/radius)/50.0)) {
        earth_x.push_back(radius*std::sin(a));
        earth_y.push_back(radius*std::cos(a));
        detector_x.push_back((radius+detector_height)*std::sin(a));
        detector_y.push_back((radius+detector_height)*std::cos(a));
        atmosphere_x.push_back((radius+atmosphere_height)*std::sin(a));
        atmosphere_y.push_back((radius+atmosphere_height)*std::cos(a));
    }
    plt.plot(earth_x,earth_y).color("b").linewidth(0.1);
    plt.plot(detector_x,detector_y).color("b").linewidth(0.1);
    plt.plot(atmosphere_x,atmosphere_y).color("b").linewidth(0.1);

    // Get the limits of the plot
    std::array<float,4> limits{min_x,x,min_y,float(radius+atmosphere_height)};
    if ((limits[1]-limits[0])*config.plot_width>(limits[3]-limits[2])*config.plot_height)
        limits[2] = limits[3] - (limits[1]-limits[0])*config.plot_height/config.plot_width;
    else
        limits[0] = limits[1] - (limits[3]-limits[2])*config.plot_width/config.plot_height;

    float d = (limits[3]-limits[2])/32.0;
    limits[0] -= d; limits[1] += d; limits[2] -=d; limits[3] += d;

    // Vertical line at x=0
    plt.plot({0,0},{0,float(radius+atmosphere_height)}).color("k").linewidth(0.5).alpha(0.5);

    // Create figure and save it.
    std::cout<<"Saving figure..."<<std::endl;
    plt.axis(limits).linewidth(0).xticks({}).yticks({})
        .figsize({config.plot_width,config.plot_height}).savefig(config.output);

}
