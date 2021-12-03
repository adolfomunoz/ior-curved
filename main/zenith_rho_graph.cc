#include <svg-cpp-plot/svg-cpp-plot.h>
#include <ivp/ivp.h>
#include <mj2/tracer/primitives/sphere.h>
#include "../eq/fermat.h"
#include <cmath>
#include <string>

#include "../atmospheres/simple.h"
#include "../atmospheres/strange.h"

/**
 * Plot [Figure 1] and [Figure 2] that shows how distance to the Z axis (rho) changes depending the zenith of the
 * Cherenkov cone's centre.
 * Also plot [Figure 8] with the flight time of the rays.
 */
int main(int argc, char** argv) {
    /** Initialization */
    /** Avaiable parameters:
     * -omega_ch <omega_ch> // Cherenkov angle
     */

    float width = 200; // Plot width
    float height = 200; // Plot height
    const double radius = 6370949; // Earth radius (m)
    const float detector_height = 0; // Meters adobe radius
    const float step_size = 100.0f; // Step size for the numeric method
    float min_zenith = 87 ;
    float max_zenith = 87.05 ;
    float zenith_step = 0.1 ;

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
    IVP::RungeKutta2 method(step_size);

    // Scene
    tracer::Sphere earth(Eigen::Vector3f(0,0,0),radius);
    tracer::Sphere surface(Eigen::Vector3f(0,0,0),radius+detector_height);




    auto time_elapsed_straight = [=] (Eigen::Vector3f ori, Eigen::Vector3f des, bool ior_at_ori = true) {
        float c = 299792458; // m/s

        Eigen::Vector3f dir = Eigen::Vector3f(des[0]-ori[0],des[1]-ori[1],des[2]-ori[2]).normalized();

        float distance = std::sqrt((ori[0]-des[0])*(ori[0]-des[0]) + (ori[1]-des[1])*(ori[1]-des[1]) + (ori[2]-des[2])*(ori[2]-des[2]));
        int n_steps = ceil(distance / step_size);
        float acc_nl = 0;
        Eigen::Vector3f pos = ori;

        // Divide straight line in steps and calculate time elapsed for each
        for (int s = 0 ; s < n_steps ; s++) {
            if (s < n_steps - 1) {
                if (ior_at_ori) {
                    acc_nl += ior(pos[0], pos[1], pos[2]) * step_size;
                }
                pos[0] = pos[0] + dir[0] * step_size;
                pos[1] = pos[1] + dir[1] * step_size;
                pos[2] = pos[2] + dir[2] * step_size;
                if (!ior_at_ori) {
                    acc_nl += ior(pos[0], pos[1], pos[2]) * step_size;
                }
            } else { // Last step
                if (ior_at_ori) {
                    acc_nl += ior(pos[0], pos[1], pos[2]) * fmod(distance,step_size);
                    break;
                }
                pos[0] = pos[0] + dir[0] * fmod(distance,step_size);
                pos[1] = pos[1] + dir[1] * fmod(distance,step_size);
                pos[2] = pos[2] + dir[2] * fmod(distance,step_size);
                acc_nl += ior(pos[0], pos[1], pos[2]) * fmod(distance, step_size);
            }
        }
        return acc_nl / c;
    };





    /** Zenith loop **/
    std::list<float> zeniths, rhos, rhos_c;
    Eigen::Array<std::list<float>,3,1> rhos_z, rhos_zc, rhos_diff, tflight, tflight_c, accs, debug_sc;
    for (float ze = min_zenith; ze < max_zenith; ze += zenith_step) {
        printf("\nZenith: %f\n", ze);
        float angle = ze * M_PI/180;
        zeniths.push_back(ze);

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
        printf("\nEXECUTION INFORMATION:\n");
        printf("Zenith: %f degrees\n", angle*(180/M_PI));
        printf("Cherenkov angle: %f degrees\n", omega_ch*(180/M_PI));
        printf("Origin of rays: %f,%f\n", x, z);
        printf("Distance Origin-North: %f m\n", distance);
        printf("Atmosphere height: %f m\n", atmosphere_height);


        /** Ray tracing */
        std::list<float> hits_x, hits_y, nohits_x, nohits_y, hits_nonlinear_x, hits_nonlinear_y;
        float zenith_s, zenith_sz, zenith_c = 0;
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
                //float rho = (*hit).point()[0] / 1000;
                rhos.push_back(rho);
                rhos_z[i].push_back(rho);

                tflight[i].push_back(time_elapsed_straight(ray.origin(), (*hit).point()));
                std::cout << "Straight ray flight time: " << tflight[i].back() << "\n";

                Eigen::Vector3f ray_dir = Eigen::Vector3f(origin[0]-(*hit).point()[0],origin[1]-(*hit).point()[1],origin[2]-(*hit).point()[2]); // Ray
                Eigen::Vector3f ref_dir = Eigen::Vector3f((*hit).point()[0],(*hit).point()[1],(*hit).point()[2]); // Reference (Earth norm)
                zenith_s = std::atan2((ray_dir.cross(ref_dir)).norm(),ray_dir.dot(ref_dir));
                std::cout << "Angle of arrival (straight ray): " << zenith_s * (180 / M_PI) << "\n";
                std::cout << "Straight HIT: " << (*hit).point()[0] << ", " << (*hit).point()[1] << ", " << (*hit).point()[2] << "\n";

                Eigen::Vector3f z_dir = Eigen::Vector3f(0,0,1); // Reference (Earth norm)
                zenith_sz = std::atan2((ray_dir.cross(z_dir)).norm(),ray_dir.dot(z_dir));
                std::cout << "Angle of shoot (straight ray): " << zenith_sz * (180 / M_PI) << "\n";

            }


            Eigen::Array<float, 6, 1> ini;
            ini(Eigen::seq(Eigen::fix<0>, Eigen::fix<2>)) = ray.origin();
            ini(Eigen::seq(Eigen::fix<3>, Eigen::fix<5>)) = ray.direction();
            Eigen::Vector3f ray1_dir = ray.direction();
            // Trace curved ray towards de Earth
            float acc_t_end = 0;
            float acc_t_start = 0;
            float tflight_sc = 0;
            for (auto s : method.steps(function, 0.0f, ini, 4.0f * distance)) {
                ray.set_range_max(s.step());
                if (auto hit = earth.trace(ray)) {
                    hits_nonlinear_x.push_back((*hit).point()[0]);
                    hits_nonlinear_y.push_back((*hit).point()[2]);
                    float rho = (std::sqrt(((*hit).point()[0] * (*hit).point()[0]))) /1000; // km
                    //float rho = (*hit).point()[0] / 1000;
                    rhos.push_back(rho);
                    rhos_zc[i].push_back(rho);

                    acc_t_start += time_elapsed_straight(ray.origin(), (*hit).point());
                    acc_t_end += time_elapsed_straight(ray.origin(), (*hit).point(), false);
                    tflight_sc = time_elapsed_straight(origin, (*hit).point());

                    Eigen::Vector3f ray_dir = Eigen::Vector3f(ray.origin()[0]-(*hit).point()[0],ray.origin()[1]-(*hit).point()[1],ray.origin()[2]-(*hit).point()[2]); // Ray
                    Eigen::Vector3f ref_dir = Eigen::Vector3f((*hit).point()[0],(*hit).point()[1],(*hit).point()[2]); // Reference (Earth norm)
                    zenith_c = std::atan2((ray_dir.cross(ref_dir)).norm(),ray_dir.dot(ref_dir));
                    std::cout << "Angle of arrival (curved ray): " << zenith_c * (180 / M_PI) << "\n";
                    std::cout << "Curved HIT: " << (*hit).point()[0] << ", " << (*hit).point()[1] << ", " << (*hit).point()[2] << "\n";
                    Eigen::Vector3f z_dir = Eigen::Vector3f(0,0,-1); // Reference Z axis
                    float zenith_cz = std::atan2((ray1_dir.cross(z_dir)).norm(),ray1_dir.dot(z_dir));
                    std::cout << "Angle of shoot (curved ray): " << zenith_cz * (180 / M_PI) << "\n";

                    // Straight ray towards hit point of the curved ray
                    Eigen::Vector3f rayX_dir = Eigen::Vector3f(origin[0]-(*hit).point()[0],origin[1]-(*hit).point()[1],origin[2]-(*hit).point()[2]); // Ray
                    float zenith_sX = std::atan2((rayX_dir.cross(ref_dir)).norm(),rayX_dir.dot(ref_dir));
                    std::cout << "Angle of arrival (straight ray to curved hit): " << zenith_sX * (180 / M_PI) << "\n";
                    z_dir = Eigen::Vector3f(0,0,1); // Reference Z axis
                    float zenith_sXz = std::atan2((rayX_dir.cross(z_dir)).norm(),rayX_dir.dot(z_dir));
                    std::cout << "Angle of shoot (straight ray to curved hit): " << zenith_sXz * (180 / M_PI) << "\n";

                    std::cout << "Difference of angle (straight - curved): " << (zenith_s-zenith_c) * (180 / M_PI) << "\n";
                    std::cout << "Difference of angle (straight_to_curved_hit - curved): " << (zenith_sX-zenith_c) * (180 / M_PI) << "\n";
                    std::cout << "Straight to curved HIT: " << (*hit).point()[0] << ", " << (*hit).point()[1] << ", " << (*hit).point()[2] << "\n";
                    std::cout << "Difference of angle (shoot) (straight - curved): " << (zenith_sz-zenith_cz) * (180 / M_PI) << "\n";
                    std::cout << "Difference of angle (shoot) (straight_to_curved_hit - curved): " << (zenith_sXz-zenith_cz) * (180 / M_PI) << "\n";

                    break;
                } else {
                    acc_t_start += time_elapsed_straight(ray.origin(), s.y()(Eigen::seq(Eigen::fix<0>, Eigen::fix<2>)));
                    acc_t_end += time_elapsed_straight(ray.origin(),  s.y()(Eigen::seq(Eigen::fix<0>, Eigen::fix<2>)), false);

                    if (std::abs(s.y()[1]) > 1.e-2)
                        std::cerr << "Warning : displacement in y :" << s.y()[1] << std::endl;

                    ray = tracer::Ray(s.y()(Eigen::seq(Eigen::fix<0>, Eigen::fix<2>)),
                                      s.y()(Eigen::seq(Eigen::fix<3>, Eigen::fix<5>)));
                }
            }

            if (isHit) {
                tflight_c[i].push_back(acc_t_end);
                std::cout << "Flight time (curved) (ior end): " << acc_t_end << "\n";
                std::cout << "Flight time (curved) (ior start): " << acc_t_start << "\n";
                std::cout << "Flight time (straight to curved hit): " << tflight_sc << "\n";
                std::cout << "Difference (straight - curved): " << tflight_sc - acc_t_start << "\n";

                if (tflight_sc - acc_t_start < 0) {
                    debug_sc[i].push_back(tflight_sc - acc_t_start);
                    std::cout << "ERROR: Straight path is faster than curved (Fermat)" << "\n";
                }
                accs[i].push_back((acc_t_end) - (acc_t_start));
                float a_ = rhos_zc[i].back();
                float b_ = rhos_z[i].back();
                float c_ = a_ - b_;
                if (i >= 2) {
                    c_ *= -1;
                }
                rhos_diff[i].push_back(c_);
                printf("[%d] Diff (%f - %f = %f)\n", i, a_, b_, c_);
            }
            std::cout << "\n";
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
    //plt.plot({zeniths.front(),zeniths.back()}, {10e-1,10e-1}).linestyle("-").color( "k").linewidth(1).alpha(0.25);
    //plt.plot({zeniths.front(),zeniths.back()}, {-10e-1,-10e-1}).linestyle("-").color( "k").linewidth(1).alpha(0.25);

    plt.xlabel("Zenith (degrees)\n");
    plt.ylabel("Distance from hit to Z axis (km)\n\n");
    plt.title("[Fig1] Red: -1 | Green: 0 | Blue: +1 | Continuous: Curved ");
    float max_y = *(std::max_element(rhos_z[2].begin(),rhos_z[2].end()));
    float min_y = *(std::min_element(rhos_zc[1].begin(),rhos_zc[1].end()));
    float min_x = *(std::min_element(zeniths.begin(), zeniths.end()));
    float max_x = *(std::max_element(zeniths.begin(), zeniths.end()));
    std::array<float,4> limits{min_x,max_x,10e-2,10e2};
    plt.set_yscale(svg_cpp_plot::symlog).base(10);
    plt.set_yticks({0,10e-3,10e-2,10e-1,10e0,10e1,10e2});//,-10e-3,-10e-2,-10e-1,-10e0,-10e1,-10e2});
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

    /**
     * Figure 8: Flight time
     */
    printf("Creating flight time figure...\n");
    svg_cpp_plot::SVGPlot plt_t;
    plt_t.plot(zeniths, tflight[0]).linestyle("--").color( "r").linewidth(1);
    plt_t.plot(zeniths, tflight[1]).linestyle("--").color( "g").linewidth(1);
    plt_t.plot(zeniths, tflight[2]).linestyle("--").color( "b").linewidth(1);
    plt_t.scatter(zeniths,tflight[0]).c("r"); // -1
    plt_t.scatter(zeniths,tflight[1]).c("g"); // 0
    plt_t.scatter(zeniths,tflight[2]).c("b"); // 1

    plt_t.plot(zeniths, tflight_c[0]).linestyle("-").color( "r").linewidth(1);
    plt_t.plot(zeniths, tflight_c[1]).linestyle("-").color( "g").linewidth(1);
    plt_t.plot(zeniths, tflight_c[2]).linestyle("-").color( "b").linewidth(1);
    plt_t.scatter(zeniths,tflight_c[0]).c("r"); // -1
    plt_t.scatter(zeniths,tflight_c[1]).c("g"); // 0
    plt_t.scatter(zeniths,tflight_c[2]).c("b"); // 1

    plt_t.plot({zeniths.front(),zeniths.back()}, {0,0}).linestyle("-").color( "k").linewidth(1);

    plt_t.xlabel("Zenith (degrees)\n");
    plt_t.ylabel("Flight time (s)\n\n");
    plt_t.title("[Fig8] Red: -1 | Green: 0 | Blue: +1 | Continuous: Curved ");
    max_y = *(std::max_element(tflight[2].begin(),tflight[2].end()));
    min_y = *(std::min_element(tflight_c[1].begin(),tflight_c[1].end()));
    min_x = *(std::min_element(zeniths.begin(), zeniths.end()));
    max_x = *(std::max_element(zeniths.begin(), zeniths.end()));
    //limits = {min_x,max_x,min_y,max_y};
    limits = {min_x,max_x,min_y,max_y};
    //plt_t.set_yscale(svg_cpp_plot::symlog).base(10);
    //plt_t.set_yticks({0,10e-3,10e-2,10e-1,10e0,10e1});
    plt_t.set_xticks({80,80.5,81,81.5,82,82.5,83,83.5,84,84.5,85,85.5,86,86.5,87,87.5,88});
    plt_t.axis(limits).linewidth(0.5).savefig("Fig8_flight_time.svg");

    /*float sum = 0;
    for (float t : debug_sc[0]) {
        sum += t;
    }
    std::cout << "Debug mean: " << sum / debug_sc[0].size();*/

    // 0.1 o 0.01
    // 8.43768333e-8 (80 ns)
    // -1.437639e-10 (-0.14 ns)

}