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
    double atmosphere_height = -10+12/cos(angle); // Km
    for (int i = 0; i<(argc-1); ++i) { // Custom angle as argument
        if (std::string(argv[i]) == "-height") atmosphere_height = atof(argv[++i]);
    }
    printf("Atmosphere height: %f km\n", atmosphere_height);
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
    float max_zenith = 90 ;
    float zenith_step = 0.17 ;
    std::list<float> zeniths, rhos, rhos_c;
    Eigen::Array<std::list<float>,3,1> rhos_z, rhos_zc, rhos_diff;
    for (float ze = min_zenith; ze < max_zenith; ze += zenith_step) {
        printf("\nZenith: %f\n", ze);
        float zenith = ze * M_PI/180;
        zeniths.push_back(ze);
        //zeniths.push_back(1/cos(zenith));

        /** Get atmosphere border (origin of the ray) from the Earth with the zenith angle */
        //To solve where the origin of the ray is (x,z coordinates) we need to solve an order 2 equation and keep the positive root

        double tana = std::tan(zenith);
        double tana2 = tana*tana;

        double a = tana2 + 1;
        double b = -2.0f*radius*tana2;
        double c = radius*radius*tana2 - (radius+atmosphere_height)*(radius+atmosphere_height);

        float z = (-b + std::sqrt(b*b - 4.0f*a*c))/(2*a);
        float x = (z-radius)*tana;

        printf("Origin of rays: %f,%f\n", x, z);
        Eigen::Vector3f origin(x,0,z);
        float distance = std::sqrt((z-radius)*(z-radius)+(x-0)*(x-0)); // distance from origin to north pole
        printf("Distance Origin-North: %f km\n", distance);

        auto function = fermat(ior,dior);
        IVP::Adaptive<IVP::Dopri> method(100,1.e-3); // (Pasos inicialmente, tolerancia)
        //IVP::Dopri method(100);
        //    IVP::Euler method(100);

        /** Cherenkov cone */
        // For -1, 0 and 1 angles
        double omega_ch = 88*M_PI/180.0; //std::acos(1/ior(x,0,z))*M_PI/180.0; //arccos(1/ior)
        for (int i = 0; i<(argc-1); ++i) { // Custom angle as argument
            if (std::string(argv[i]) == "-omega_ch") omega_ch = atof(argv[++i])*M_PI/180.0;
        }
        printf("omega_ch: %f degrees\n", omega_ch*(180/M_PI));

        int i = 0;
        for (float angle = (zenith-omega_ch); angle<=(zenith+1.5*omega_ch); angle+=omega_ch) {
            bool isHit = false;
            printf("Angle: %f\n", angle);
            std::list<float> hits_x, hits_y, hits_nonlinear_x, hits_nonlinear_y;
            tracer::Ray ray(origin, Eigen::Vector3f(-std::sin(angle), 0, -std::cos(angle)));
            // Trace ray towards the Earth
            if (auto hit = earth.trace(ray)) {
                isHit = true;
                hits_x.push_back((*hit).point()[0]);
                hits_y.push_back((*hit).point()[2]);
                float rho = std::sqrt(((*hit).point()[0] * (*hit).point()[0]));
                rhos.push_back(rho);
                rhos_z[i].push_back(rho);
                //printf("Hit (rho: %f)\n",rho);
                //printf("Hit: %f,%f (rho: %f)\n", (*hit).point()[0], (*hit).point()[2], rho);
            }
            Eigen::Array<float, 6, 1> ini;
            ini(Eigen::seq(Eigen::fix<0>, Eigen::fix<2>)) = ray.origin();
            ini(Eigen::seq(Eigen::fix<3>, Eigen::fix<5>)) = ray.direction();
            // Trace curved ray towards de Earth
            for (auto s : method.steps(function, 0.0f, ini, 1.5f * distance)) {//1.0f*float(atmosphere_height))) {
                ray.set_range_max(s.step());
                if (auto hit = earth.trace(ray)) {
                    isHit = true;
                    hits_nonlinear_x.push_back((*hit).point()[0]);
                    hits_nonlinear_y.push_back((*hit).point()[2]);
                    float rho = std::sqrt(((*hit).point()[0] * (*hit).point()[0]));
                    rhos_c.push_back(rho);
                    rhos_zc[i].push_back(rho);
                    //printf("Hit Curved (rho: %f)\n",rho);
                    //printf("Hit Curved: %f,%f (rho: %f)\n", (*hit).point()[0], (*hit).point()[2], rho);
                    break;
                } else {
                    if (std::abs(s.y()[1]) > 1.e-2)
                        std::cerr << "Warning : displacement in y :" << s.y()[1] << std::endl;
                    ray = tracer::Ray(s.y()(Eigen::seq(Eigen::fix<0>, Eigen::fix<2>)),
                                      s.y()(Eigen::seq(Eigen::fix<3>, Eigen::fix<5>)));
                }
            }
            if (isHit) {
                rhos_diff[i].push_back(rhos_zc[i].back() - rhos_z[i].back());
            } else {
                rhos_diff[i].push_back(-10000);
            }

            i++;
        }
    }

    printf("Creating figure...");
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

    plt.xlabel("Zenith (degrees)\n");
    plt.ylabel("Distance from hit to Z axis\n\n");
    plt.title("Red: -1 | Green: 0 | Blue: +1 | Continuous: Curved ");
    float max_y = *(std::max_element(rhos_zc[0].begin(),rhos_zc[0].end()));
    std::array<float,4> limits{79.5,90.5,0,max_y+100};
    plt.axis(limits).linewidth(0.5).savefig("rho_graph.svg");

    printf("Creating diff figure...");
    svg_cpp_plot::SVGPlot plt_diff;
    plt_diff.plot(zeniths, rhos_diff[0]).linestyle("-").color( "r").linewidth(1);
    plt_diff.plot(zeniths, rhos_diff[1]).linestyle("-").color( "g").linewidth(1);
    plt_diff.plot(zeniths, rhos_diff[2]).linestyle("-").color( "b").linewidth(1);
    plt_diff.scatter(zeniths,rhos_diff[0]).c("r"); // -1
    plt_diff.scatter(zeniths,rhos_diff[1]).c("g"); // 0
    plt_diff.scatter(zeniths,rhos_diff[2]).c("b"); // 1

    plt_diff.xlabel("Zenith (degrees)\n");
    plt_diff.ylabel("Difference of rho between straight and curved");
    plt_diff.title("Red: -1 | Green: 0 | Blue: +1 ");
    max_y = *(std::max_element(rhos_diff[1].begin(),rhos_diff[1].end()));
    limits = {79.5,90.5,0,max_y+100};
    plt_diff.axis(limits).linewidth(0.5).savefig("rho_graph_diff.svg");

}