#include <svg-cpp-plot/svg-cpp-plot.h>
#include <ivp/ivp.h>
#include <mj2/tracer/primitives/sphere.h>
#include "../eq/fermat.h"
#include <cmath>
#include <string>
/**
 * Plot [Figure 1] and [Figure 2] that shows how distance to the Z axis (rho) changes depending the zenith of the
 * Cherenkov cone's centre.
 */
int main(int argc, char** argv) {
    /** Parameters */
    // Plot
    float width = 200;
    float height = 50;

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

    float min_zenith = 80 ;
    float max_zenith = 88.05 ;
    float zenith_step = 0.1 ;
    std::list<float> zeniths, rhos, rhos_c;
    Eigen::Array<std::list<float>,3,1> rhos_z, rhos_zc, rhos_diff;
    for (float ze = min_zenith; ze < max_zenith; ze += zenith_step) {
        printf("\nZenith: %f\n", ze);
        float angle = ze * M_PI/180;
        zeniths.push_back(ze);

        /** Get atmosphere border (origin of the ray) from the Earth with the zenith angle */
        //To solve where the origin of the ray is (x,z coordinates) we need to solve an order 2 equation and keep the positive root
        float distance = (-10 + 12 / std::cos(angle)) * 1000; // x1000 km -> m
        float x = distance * std::sin(angle);
        float z = (distance * std::cos(angle)) + radius;
        float atmosphere_height = std::sqrt(x * x + z * z) - radius;

        printf("Origin of rays: %f,%f\n", x, z);
        Eigen::Vector3f origin(x, 0, z);
        printf("Distance Origin-North: %f m\n", distance);
        printf("Atmosphere height: %f m\n", atmosphere_height);

        auto function = fermat(ior, dior);
        //IVP::Adaptive<IVP::Dopri> method(100,1.e-3); // (Pasos inicialmente, tolerancia)
        IVP::RungeKutta2 method(2000.0f);

        /** Cherenkov cone */
        // For -1, 0 and 1 angles
        double omega_ch = 88 * M_PI / 180.0; //std::acos(1/ior(x,0,z))*M_PI/180.0; //arccos(1/ior)
        for (int i = 0; i < (argc - 1); ++i) { // Custom angle as argument
            if (std::string(argv[i]) == "-omega_ch") omega_ch = atof(argv[++i]) * M_PI / 180.0;
        }
        printf("omega_ch: %f degrees\n", omega_ch * (180 / M_PI));

        /** Plotting */
        svg_cpp_plot::SVGPlot plt;
        std::list<float> hits_x, hits_y, nohits_x, nohits_y, hits_nonlinear_x, hits_nonlinear_y;
        int i = 0;
        for (float a = (angle - omega_ch); a <= (angle + 1.5 * omega_ch); a += omega_ch) {
            bool isHit = false;
            tracer::Ray ray(origin, Eigen::Vector3f(-std::sin(a), 0, -std::cos(a)));
            // Trace ray towards the Earth
            if (auto hit = earth.trace(ray)) {
                isHit = true;
                hits_x.push_back((*hit).point()[0]);
                hits_y.push_back((*hit).point()[2]);
                float rho = (std::sqrt(((*hit).point()[0] * (*hit).point()[0])))/1000; // km
                rhos.push_back(rho);
                rhos_z[i].push_back(rho);
            }
            Eigen::Array<float, 6, 1> ini;
            ini(Eigen::seq(Eigen::fix<0>, Eigen::fix<2>)) = ray.origin();
            ini(Eigen::seq(Eigen::fix<3>, Eigen::fix<5>)) = ray.direction();
            // Trace curved ray towards de Earth
            for (auto s : method.steps(function, 0.0f, ini, 4.0f * distance)) {//1.0f*float(atmosphere_height))) {
                ray.set_range_max(s.step());
                if (auto hit = earth.trace(ray)) {
                    hits_nonlinear_x.push_back((*hit).point()[0]);
                    hits_nonlinear_y.push_back((*hit).point()[2]);
                    float rho = (std::sqrt(((*hit).point()[0] * (*hit).point()[0]))) /1000; // km
                    rhos.push_back(rho);
                    rhos_zc[i].push_back(rho);
                    break;
                } else {
                    if (std::abs(s.y()[1]) > 1.e-2)
                        std::cerr << "Warning : displacement in y :" << s.y()[1] << std::endl;
                    ray = tracer::Ray(s.y()(Eigen::seq(Eigen::fix<0>, Eigen::fix<2>)),
                                      s.y()(Eigen::seq(Eigen::fix<3>, Eigen::fix<5>)));
                }
            }
            if (isHit) {
                float a = rhos_zc[i].back();
                float b = rhos_z[i].back();
                float c = a - b;
                if (i >= 2) {
                    c *= -1;
                }
                rhos_diff[i].push_back(c);
                printf("[%d] Diff (%f - %f = %f)\n", i, a, b, c);
            }

            i++;
        }
    }

    /**
     * Figure 1: Distance from hit to Z axis
     */
    printf("Creating figure...\n");
    svg_cpp_plot::SVGPlot plt;
    plt.plot(zeniths, rhos_z[0]).linestyle("--").color( "r").linewidth(1);
    plt.plot(zeniths, rhos_z[1]).linestyle("--").color( "g").linewidth(1);
    plt.plot(zeniths, rhos_z[2]).linestyle("--").color( "b").linewidth(1);
    plt.scatter(zeniths,rhos_z[0]).c("r"); // -1
    plt.scatter(zeniths,rhos_z[1]).c("g"); // 0
    plt.scatter(zeniths,rhos_z[2]).c("b"); // 1

    plt.plot(zeniths, rhos_zc[0]).linestyle("-").color( "r").linewidth(1);
    plt.plot(zeniths, rhos_zc[1]).linestyle("-").color( "g").linewidth(1);
    plt.plot(zeniths, rhos_zc[2]).linestyle("-").color( "b").linewidth(1);
    plt.scatter(zeniths,rhos_zc[0]).c("r"); // -1
    plt.scatter(zeniths,rhos_zc[1]).c("g"); // 0
    plt.scatter(zeniths,rhos_zc[2]).c("b"); // 1

    plt.plot({zeniths.front(),zeniths.back()}, {0,0}).linestyle("-").color( "k").linewidth(1);

    plt.xlabel("Zenith (degrees)\n");
    plt.ylabel("Distance from hit to Z axis (km)\n\n");
    plt.title("[Fig1] Red: -1 | Green: 0 | Blue: +1 | Continuous: Curved ");
    float max_y = *(std::max_element(rhos_z[2].begin(),rhos_z[2].end()));
    float min_y = *(std::min_element(rhos_zc[1].begin(),rhos_zc[1].end()));
    float min_x = *(std::min_element(zeniths.begin(), zeniths.end()));
    float max_x = *(std::max_element(zeniths.begin(), zeniths.end()));
    std::array<float,4> limits{min_x,max_x,0.1,10e2};
    plt.set_yscale(svg_cpp_plot::symlog).base(10);
    plt.set_yticks({0,10e-3,10e-2,10e-1,10e0,10e1,10e2});
    plt.set_xticks({80,80.5,81,81.5,82,82.5,83,83.5,84,84.5,85,85.5,86,86.5,87,87.5,88});
    plt.axis(limits).linewidth(0.5).savefig("Fig1_rho_graph.svg");

    /**
     * Figure 2: Difference between straight and curved rays in the distance from hit to the Z axis.
     */
    printf("Creating diff figure...\n");
    svg_cpp_plot::SVGPlot plt_diff;
    plt_diff.plot(zeniths, rhos_diff[0]).linestyle("-").color( "r").linewidth(1);
    plt_diff.plot(zeniths, rhos_diff[1]).linestyle("-").color( "g").linewidth(1);
    plt_diff.plot(zeniths, rhos_diff[2]).linestyle("-").color( "b").linewidth(1);

    plt_diff.scatter(zeniths,rhos_diff[0]).c("r"); // -1
    plt_diff.scatter(zeniths,rhos_diff[1]).c("g"); // 0
    plt_diff.scatter(zeniths,rhos_diff[2]).c("b"); // 1
    plt_diff.plot({zeniths.front(),zeniths.back()}, {0,0}).linestyle("-").color( "k").linewidth(1);

    plt_diff.xlabel("Zenith (degrees)\n");
    plt_diff.ylabel("Difference of rho between straight and curved (km)");
    plt_diff.title("[Fig2] Red: -1 | Green: 0 | Blue: +1 ");
    max_y = *(std::max_element(rhos_diff[1].begin(),rhos_diff[1].end()));
    min_y = *(std::min_element(rhos_diff[2].begin(),rhos_diff[2].end()));
    min_x = *(std::min_element(zeniths.begin(), zeniths.end()));
    max_x = *(std::max_element(zeniths.begin(), zeniths.end()));
    //limits = {min_x,max_x,min_y,max_y};
    limits = {80,88,0,10e1};
    plt_diff.set_yscale(svg_cpp_plot::symlog).base(10);
    plt_diff.set_yticks({0,10e-3,10e-2,10e-1,10e0,10e1});
    plt_diff.set_xticks({80,80.5,81,81.5,82,82.5,83,83.5,84,84.5,85,85.5,86,86.5,87,87.5,88});
    plt_diff.axis(limits).linewidth(0.5).savefig("Fig2_rho_graph_diff.svg");
}