#include <svg-cpp-plot/svg-cpp-plot.h>
#include <ivp/ivp.h>
#include <mj2/tracer/primitives/sphere.h>
#include "../eq/fermat.h"
#include <cmath>
#include <string>
/**
 * Plot [Figure 3], [Figure 4], [Figure 5] and [Figure 6] using rho (distance to the Z axis) of the ray hits inside the
 * Cherenkov cone.
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

    /** Get atmosphere border (origin of the ray) from the Earth with the zenith angle */
    //To solve where the origin of the ray is (x,z coordinates) we need to solve an order 2 equation and keep the positive root
    float distance = (-10 + 12 / std::cos(angle)) * 1000;
    float x = distance * std::sin(angle);
    float z = ( distance * std::cos(angle) ) + radius;
    float atmosphere_height = std::sqrt(x*x + z*z) - radius;

    printf("Origin of rays: %f,%f\n", x, z);
    Eigen::Vector3f origin(x,0,z);
    printf("Distance Origin-North: %f km\n", distance);
    printf("Atmosphere height: %f km\n", atmosphere_height);

    auto function = fermat(ior,dior);
    //IVP::Adaptive<IVP::Dopri> method(100,1.e-3); // (Initial steps, tolerance)
    IVP::RungeKutta2 method(2000.0f);

    /** Cherenkov cone */
    // For -1, 0 and 1 angles
    double omega_ch = 88*M_PI/180.0;
    for (int i = 0; i<(argc-1); ++i) { // Custom angle as argument
        if (std::string(argv[i]) == "-omega_ch") omega_ch = atof(argv[++i])*M_PI/180.0;
    }
    printf("omega_ch: %f degrees\n", omega_ch*(180/M_PI));

    std::list<float> zeniths, angles, omegas, omegas_c, rhos, rhos_c, rhos_diff;
    float angle_step = 0.03 * M_PI/180.0;
    int i = 0;
    for (float a = (angle-omega_ch); a<=(angle+omega_ch+0.9*angle_step); a+=angle_step) {
        printf("Angle: %f\n", a * (180/M_PI));
        zeniths.push_back(a *(180/M_PI));
        angles.push_back(((angle*(180/M_PI)) - (a *(180/M_PI))) / (omega_ch*(180/M_PI)));
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
        for (auto s : method.steps(function, 0.0f, ini, 4.0f * distance)) {//1.0f*float(atmosphere_height))) {
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
                printf("Curved hit: rho_c: %f | omega_c: %f\n", rhos_c.back(), omega*(180/M_PI));
                break;
            } else {
                if (std::abs(s.y()[1]) > 1.e-2)
                    std::cerr << "Warning : displacement in y :" << s.y()[1] << std::endl;
                ray = tracer::Ray(s.y()(Eigen::seq(Eigen::fix<0>, Eigen::fix<2>)),
                                  s.y()(Eigen::seq(Eigen::fix<3>, Eigen::fix<5>)));
            }
        }

        rhos_diff.push_back(rhos_c.back() - rhos.back());
        i++;
    }

    /**
    * Figure 3: Graph with the angle to the Z axis (X axis) and the distance to the Z axis (Y axis)
    */
    printf("Creating Zenith figure...\n");
    svg_cpp_plot::SVGPlot plt_zenith;
    plt_zenith.plot(zeniths, rhos).linestyle("-").color( "r").linewidth(1);
    plt_zenith.scatter(zeniths,rhos).c("r").s(2).alpha(0.5);
    plt_zenith.plot(zeniths, rhos_c).linestyle("-").color( "b").linewidth(1);
    plt_zenith.scatter(zeniths,rhos_c).c("b").s(2).alpha(0.5);
    plt_zenith.plot({zeniths.front(),zeniths.back()}, {0,0}).linestyle("-").color( "k").linewidth(1);
    printf("debug %f", rhos_c.back());

    plt_zenith.xlabel("Zenith (degrees)\n");
    plt_zenith.ylabel("Distance from hit to Z axis\n\n");
    plt_zenith.title("Red: Hits of straight rays | Blue: Hits of curved rays");
    float min_x = *(std::min_element(zeniths.begin(),zeniths.end()));
    float min_y = *(std::min_element(rhos.begin(),rhos.end()));
    float max_x = *(std::max_element(zeniths.begin(),zeniths.end()));
    float max_y = *(std::max_element(rhos.begin(),rhos.end()));
    std::array<float,4> limits{min_x,max_x,min_y,max_y};
    plt_zenith.axis(limits).linewidth(1).savefig("Fig3_rho_ch_cone_zenith.svg");

    /**
     * Figure 4: Graph with the ray (-1, 0, 1) (X axis) and angle to the normal (Y axis)
     */
    printf("Creating Normal figure...\n");
    svg_cpp_plot::SVGPlot plt_normal;
    plt_normal.plot(angles, omegas).linestyle("-").color( "r").linewidth(1);
    plt_normal.scatter(angles,omegas).c("r").s(2).alpha(0.5);
    plt_normal.plot(angles, omegas_c).linestyle("-").color( "b").linewidth(1);
    plt_normal.scatter(angles,omegas_c).c("b").s(2).alpha(0.5);
    plt_normal.plot({angles.front(),angles.back()}, {0,0}).linestyle("-").color( "k").linewidth(1);

    plt_normal.xlabel("Cherenkov ray\n");
    plt_normal.ylabel("Angle to the normal (degrees)\n\n");
    plt_normal.title("[Fig4] Red: Hits of straight rays | Blue: Hits of curved rays\n");
    min_x = *(std::min_element(angles.begin(),angles.end()));
    min_y = *(std::min_element(omegas.begin(),omegas.end()));
    max_x = *(std::max_element(angles.begin(),angles.end()));
    max_y = *(std::max_element(omegas.begin(),omegas.end()));
    limits = {min_x,max_x,min_y,max_y};
    plt_normal.axis(limits).linewidth(1).savefig("Fig4_rho_ch_cone_normal.svg");

    /**
     * Figure 5: Graph with angle to the Z axis (X axis) and difference between rays of the distance to the
     * Z axis (rho) (Y axis)
     */
    printf("Creating Zenith Diff figure...\n");
    svg_cpp_plot::SVGPlot plt_zenith_d;
    plt_zenith_d.plot(zeniths, rhos_diff).linestyle("-").color( "r").linewidth(1);
    plt_zenith_d.scatter(zeniths,rhos_diff).c("r").s(2).alpha(0.5);
    plt_zenith_d.plot({zeniths.front(),zeniths.back()}, {0,0}).linestyle("-").color( "k").linewidth(1);

    plt_zenith_d.xlabel("Zenith (degrees)\n");
    plt_zenith_d.ylabel("Difference of rho between straight and curved");
    min_x = *(std::min_element(zeniths.begin(),zeniths.end()));
    min_y = *(std::min_element(rhos_diff.begin(),rhos_diff.end()));
    max_x = *(std::max_element(zeniths.begin(),zeniths.end()));
    max_y = *(std::max_element(rhos_diff.begin(),rhos_diff.end()));
    limits = {min_x,max_x,min_y,max_y};
    plt_zenith_d.axis(limits).linewidth(1).savefig("Fig5_rho_ch_cone_zenith_diff.svg");

    /**
     * Figure 6: Graph with angle to the normal (X axis) and distance to the Z axis (rho) (Y axis)
     */
    printf("Creating Normal 2 figure...\n");
    svg_cpp_plot::SVGPlot plt_normal_2;
    plt_normal_2.plot(omegas, rhos).linestyle("-").color( "r").linewidth(1);
    plt_normal_2.scatter(omegas,rhos).c("r").s(2).alpha(0.5);
    plt_normal_2.plot(omegas_c, rhos_c).linestyle("-").color( "b").linewidth(1);
    plt_normal_2.scatter(omegas_c,rhos_c).c("b").s(2).alpha(0.5);
    plt_normal_2.plot({omegas.front(),omegas.back()}, {0,0}).linestyle("-").color( "k").linewidth(1);

    plt_normal_2.xlabel("Angle to the normal (degrees)\n");
    plt_normal_2.ylabel("Distance to the Z axis");
    min_x = *(std::min_element(omegas.begin(),omegas.end()));
    min_y = *(std::min_element(rhos.begin(),rhos.end()));
    max_x = *(std::max_element(omegas.begin(),omegas.end()));
    max_y = *(std::max_element(rhos.begin(),rhos.end()));
    limits = {min_x,max_x,min_y,max_y};
    plt_normal_2.axis(limits).linewidth(1).savefig("Fig6_rho_ch_cone_normal_2.svg");


}