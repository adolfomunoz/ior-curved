#include <svg-cpp-plot/svg-cpp-plot.h>
#include <ivp/ivp.h>
#include <mj2/tracer/primitives/sphere.h>
#include "../eq/fermat.h"
#include <cmath>
#include <string>

#include "../atmospheres/simple.h"
#include "../atmospheres/strange.h"
#include "../atmospheres/gutierrez2006.h"

/**
 * Plot of the IOR of an atmosphere.
 */
int main(int argc, char** argv) {
    auto atmospheres = pattern::make_from_commandline<std::list<Atmosphere>>(argc,argv);
    float width = 400; // Plot width
    float height = 400; // Plot height
    const float radius = 6370949; // Earth radius (m)
    const float atmosphere_height = 34000; // Atmosphere height for plotting
    const char* output = "ior.svg";
    
    std::cerr<<"Number of atmospheres = "<<atmospheres.size()<<std::endl;
 
    /** Plotting ior **/
    svg_cpp_plot::SVGPlot plt;
    for (auto atmosphere : atmospheres) {
        std::list<float> iors, heights;
        for (float h = 0; h<atmosphere_height; h+=(atmosphere_height/height)) {
            heights.push_back(h);
            iors.push_back(atmosphere.ior(0,0,radius+h));
        }
        plt.plot(iors,heights);
    }

    if (!atmospheres.empty())
        plt.axis({1,1.0005,0,atmosphere_height}).xlabel("IOR").ylabel("Height (m)").xticks({1,1.0005}).figsize({width,height}).savefig(output);
}
