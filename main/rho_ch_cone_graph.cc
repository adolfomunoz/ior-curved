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

    /** Get atmosphere border (origin of the ray) from the Earth with the zenith angle */
    //To solve where the origin of the ray is (x,z coordinates) we need to solve an order 2 equation and keep the positive root

    double tana = std::tan(angle);
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
    IVP::Adaptive<IVP::Dopri> method(100,1.e-3); // (Initial steps, tolerance)
    //IVP::Dopri method(100);
    //    IVP::Euler method(100);

    /** Cherenkov cone */
    // For -1, 0 and 1 angles
    double omega_ch = 88*M_PI/180.0;
    for (int i = 0; i<(argc-1); ++i) { // Custom angle as argument
        if (std::string(argv[i]) == "-omega_ch") omega_ch = atof(argv[++i])*M_PI/180.0;
    }
    printf("omega_ch: %f degrees\n", omega_ch*(180/M_PI));

    std::list<float> zeniths, omegas, omegas_c, rhos, rhos_c;
    Eigen::Array<std::list<float>,3,1> rhos_z, rhos_zc;
    float angle_step = 0.05 * M_PI/180.0;
    int i = 0;
    for (float a = (angle-omega_ch); a<=(angle+1.1*omega_ch); a+=angle_step) {
        zeniths.push_back(a *(180/M_PI));
        std::list<float> hits_x, hits_y, hits_nonlinear_x, hits_nonlinear_y;
        tracer::Ray ray(origin, Eigen::Vector3f(-std::sin(a), 0, -std::cos(a)));
        // Trace ray towards the Earth
        if (auto hit = earth.trace(ray)) {
            hits_x.push_back((*hit).point()[0]);
            hits_y.push_back((*hit).point()[2]);

            float rho = (*hit).point()[0];
            rhos.push_back(rho);

            float xa = (*hit).point()[0];
            float ya = (*hit).point()[2];
            float xb = x - (*hit).point()[0];
            float yb = z - (*hit).point()[2];
            float omega = std::acos((xa * xb + ya * yb) / (std::sqrt(xa*xa + ya*ya) * std::sqrt(xb*xb + yb*yb)));
            omegas.push_back(omega*(180/M_PI));
        }
        Eigen::Array<float, 6, 1> ini;
        ini(Eigen::seq(Eigen::fix<0>, Eigen::fix<2>)) = ray.origin();
        ini(Eigen::seq(Eigen::fix<3>, Eigen::fix<5>)) = ray.direction();
        // Trace curved ray towards de Earth
        for (auto s : method.steps(function, 0.0f, ini, 1.5f * distance)) {//1.0f*float(atmosphere_height))) {
            ray.set_range_max(s.step());
            if (auto hit = earth.trace(ray)) {
                hits_nonlinear_x.push_back((*hit).point()[0]);
                hits_nonlinear_y.push_back((*hit).point()[2]);

                float rho = (*hit).point()[0];
                rhos_c.push_back(rho);

                float xa = (*hit).point()[0];
                float ya = (*hit).point()[2];
                float xb = x - (*hit).point()[0];
                float yb = z - (*hit).point()[2];
                float omega = std::acos((xa * xb + ya * yb) / (std::sqrt(xa*xa + ya*ya) * std::sqrt(xb*xb + yb*yb)));
                omegas_c.push_back(omega*(180/M_PI));
                break;
            } else {
                if (std::abs(s.y()[1]) > 1.e-2)
                    std::cerr << "Warning : displacement in y :" << s.y()[1] << std::endl;
                ray = tracer::Ray(s.y()(Eigen::seq(Eigen::fix<0>, Eigen::fix<2>)),
                                  s.y()(Eigen::seq(Eigen::fix<3>, Eigen::fix<5>)));
            }
        }
        i++;
    }

    printf("Creating Zenith figure...\n");
    svg_cpp_plot::SVGPlot plt_zenith;
    plt_zenith.plot(zeniths, rhos).linestyle("-").color( "r").linewidth(1);
    plt_zenith.scatter(zeniths,rhos).c("r").s(2).alpha(0.5);
    plt_zenith.plot(zeniths, rhos_c).linestyle("-").color( "b").linewidth(1);
    plt_zenith.scatter(zeniths,rhos_c).c("b").s(2).alpha(0.5);
    plt_zenith.plot({zeniths.front(),zeniths.back()}, {0,0}).linestyle("-").color( "k").linewidth(1);

    plt_zenith.xlabel("Zenith (degrees)\n");
    plt_zenith.ylabel("Distance from hit to Z axis\n\n");
    plt_zenith.title("Red: Hits of straight rays | Blue: Hits of curved rays");
    plt_zenith.linewidth(1).savefig("rho_ch_cone_zenith.svg");


    printf("Creating Normal figure...\n");
    svg_cpp_plot::SVGPlot plt_normal;
    plt_normal.plot(omegas, rhos).linestyle("-").color( "r").linewidth(1);
    plt_normal.scatter(omegas,rhos).c("r").s(2).alpha(0.5);
    plt_normal.plot(omegas_c, rhos_c).linestyle("-").color( "b").linewidth(1);
    plt_normal.scatter(omegas_c,rhos_c).c("b").s(2).alpha(0.5);
    plt_normal.plot({zeniths.front(),zeniths.back()}, {0,0}).linestyle("-").color( "k").linewidth(1);

    plt_normal.xlabel("Angle to the normal (degrees)\n");
    plt_normal.ylabel("Distance from hit to Z axis\n\n");
    plt_normal.title("Red: Hits of straight rays | Blue: Hits of curved rays");
    plt_normal.linewidth(1).savefig("rho_ch_cone_normal.svg");


}