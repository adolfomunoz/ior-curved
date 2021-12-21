#include <svg-cpp-plot/svg-cpp-plot.h>
#include <ivp/ivp.h>
#include <mj2/tracer/primitives/sphere.h>
#include "../eq/fermat.h"
#include <cmath>
#include <string>

#include "../atmospheres/simple.h"
#include "../atmospheres/strange.h"
#include "../atmospheres/gutierrez2006.h"


struct Config : public pattern::Reflectable<Config> {
    float plot_width = 400;
    float plot_height = 400;
    float atmosphere_height = 34000; 
    std::string output = "ior.svg";
    
    auto reflect() { return std::tie(plot_width,plot_height,atmosphere_height,output); }
    auto reflect_names() const { return std::tuple("plot-width","plot_height","max-height","output"); }
    static const char* type_name() { return "config"; } 
};
/**
 * Plot of the IOR of an atmosphere.
 */
int main(int argc, char** argv) {
    auto atmospheres = pattern::make_from_commandline<std::list<Atmosphere>>(argc,argv);
    const float radius = 6370949; // Earth radius (m)
    
    Config config;
    pattern::load_commandline(config,argc,argv);

    
    std::cerr<<"Number of atmospheres = "<<atmospheres.size()<<std::endl;
 
    /** Plotting ior **/
    svg_cpp_plot::SVGPlot plt;
    for (auto atmosphere : atmospheres) {
        std::list<float> iors, heights;
        for (float h = 0; h<config.atmosphere_height; h+=(config.atmosphere_height/config.plot_height)) {
            heights.push_back(h);
            iors.push_back(atmosphere.ior(0,0,radius+h));
        }
        plt.plot(iors,heights);
    }

    if (!atmospheres.empty())
        plt.axis({1,1.0005,0,config.atmosphere_height}).xlabel("IOR").ylabel("Height (m)").xticks({1,1.0005}).figsize({config.plot_width,config.plot_height}).savefig(config.output);
}
