#include <svg-cpp-plot/svg-cpp-plot.h>
#include <ivp/ivp.h>
#include <mj2/tracer/primitives/sphere.h>
#include "../eq/fermat.h"
#include <cmath>
#include <string>

#include "../atmospheres/simple.h"
#include "../atmospheres/strange.h"

/**
 * Plot [Figure 7] that shows in overhead view (XY axis) the hits of the Cherenkov cone using the three dimensions.
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
    const float method_steps = 1000.0f; // Steps for the numeric method
    auto zeniths = {87.6}; //{82,83,84,85,86};
    auto alphas = Eigen::Vector<float,5>(1.0,0.8,0.6,0.4,0.2); // Alpha value for each zenith plotted

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
    // [ Inside the loop ]

    // Cherenkov cone. For -1, 0 and 1 angles
    float omega_ch = 1*M_PI/180.0; // Default 1º
    for (int i = 0; i<(argc-1); ++i) { // Custom angle as argument
        if (std::string(argv[i]) == "-omega_ch") omega_ch = atof(argv[++i])*M_PI/180.0;
    }

    // Print initial information
    printf("Cherenkov angle: %f degrees\n", omega_ch*(180/M_PI));






    /** Ray tracing (& path plotting) **/

    svg_cpp_plot::SVGPlot plt_xy;
    float g_max_x, g_max_y, g_min_x, g_min_y;
    int alpha = 0;
    std::list<Eigen::Array<float,3,1>> hits_curved, hits_straight;
    for (float zenith : zeniths) {
        printf("Zenith: %f\n", zenith);
        angle = zenith * (M_PI / 180.0);

        // Get atmosphere border (origin of the ray) from the Earth with the zenith angle
        float distance = (-23.3f + 19.45f/std::cos(0.99273f*angle)) * 1000; // x1000 km -> m
        //float distance = 279.86711f * 1000; // Fixed value for testing
        float x = distance * std::sin(angle);
        float y = 0;
        float z = ( distance * std::cos(angle) ) + radius + detector_height;
        float atmosphere_height = std::sqrt(x*x + y*y + z*z) - radius;
        Eigen::Vector3f origin(x,0,z);

        printf("Origin of rays: %f,%f\n", x, z);
        printf("Distance Origin-North: %f m\n", distance);
        printf("Atmosphere height: %f m\n", atmosphere_height);


        // Station basis
        Eigen::Vector3f o = Eigen::Vector3f(origin[0], origin[1], origin[2]);
        Eigen::Vector3f i = Eigen::Vector3f(-std::sin(angle), 0, -std::cos(angle));
        i = i.normalized();
        Eigen::Vector3f j = Eigen::Vector3f(1, 1, (i[0] + i[1]) / -i[2]);
        j = j.normalized();
        Eigen::Vector3f k = i.cross(j);
        k = k.normalized();

        std::cout << i.dot(j) << "\n";
        std::cout << i.dot(k) << "\n";
        std::cout << j.dot(k) << "\n";

        //printf("i: %f, %f, %f\n", i[0], i[1], i[2]);
        //printf("j: %f, %f, %f\n", j[0], j[1], j[2]);
        //printf("k: %f, %f, %f\n", k[0], k[1], k[2]);
        // Change of basis matrix
        Eigen::Matrix3f m;
        m << i[0], j[0], k[0],
                i[1], j[1], k[1],
                i[2], j[2], k[2];

        // Center ray towards Earth
        Eigen::Vector3f dir_c = m * Eigen::Vector3f(1, 0, 0);
        printf("Center ray direction: %f, %f, %f\n", dir_c[0], dir_c[1], dir_c[2]);
        tracer::Ray center_ray(origin, dir_c);
        Eigen::Vector3f center_hit;

        if (auto hit = earth.trace(center_ray)) {
            center_hit = Eigen::Vector3f((*hit).point()[0], (*hit).point()[1], (*hit).point()[2]);
            printf("Center Hit: X:%f Y:%f Z:%f\n", (*hit).point()[0], (*hit).point()[1], (*hit).point()[2]);
            printf("Center Hit (z=z-radius): X:%f Y:%f Z:%f\n", (*hit).point()[0], (*hit).point()[1],
                   (*hit).point()[2] - radius);
        }
        Eigen::Array<float, 6, 1> ini;
        ini(Eigen::seq(Eigen::fix<0>, Eigen::fix<2>)) = center_ray.origin();
        ini(Eigen::seq(Eigen::fix<3>, Eigen::fix<5>)) = center_ray.direction();
        Eigen::Vector3f center_hit_c;
        // Trace curved ray towards de Earth
        for (auto s : method.steps(function, 0.0f, ini, 4.0f * distance)) {//1.0f*float(atmosphere_height))) {
            center_ray.set_range_max(s.step());
            if (auto hit = earth.trace(center_ray)) {
                center_hit_c = Eigen::Vector3f((*hit).point()[0], (*hit).point()[1], (*hit).point()[2]);
                printf("Center curved hit: X:%f Y:%f Z:%f (Zr:%f)\n", (*hit).point()[0], (*hit).point()[1],
                       (*hit).point()[2], (*hit).point()[2] - radius);
                break;
            } else {
                center_ray = tracer::Ray(s.y()(Eigen::seq(Eigen::fix<0>, Eigen::fix<2>)),
                                         s.y()(Eigen::seq(Eigen::fix<3>, Eigen::fix<5>)));
            }
        }

        std::list<float> hits_x, hits_y, hits_z, hits_nonlinear_x, hits_nonlinear_y, hits_nonlinear_z, nohits_x, nohits_y, path_x, path_y, path_z, nohits_z;
        std::list<float> xc, yc, zc;
        for (float c = 0; c < 2 * M_PI; c += M_PI / 32) {
            // Cone calculation
            Eigen::Vector3f v = Eigen::Vector3f(std::cos(omega_ch), std::sin(omega_ch) * std::cos(c),
                                                std::sin(omega_ch) * std::sin(c));
            Eigen::Vector3f dir = m * v;
            xc.push_back(dir[0]);
            yc.push_back(dir[1]);
            zc.push_back(dir[2]);

            // Trace ray towards the Earth
            tracer::Ray ray(origin, dir);
            if (auto hit = earth.trace(ray)) {
                hits_straight.push_back( Eigen::Array<float,3,1>((*hit).point()[0], (*hit).point()[1], (*hit).point()[2]));
                hits_x.push_back((*hit).point()[0]);
                hits_y.push_back((*hit).point()[1]);
                hits_z.push_back((*hit).point()[2]);
                printf("Hit: X:%f Y:%f Z:%f (Zr:%f) (%f)\n", (*hit).point()[0], (*hit).point()[1], (*hit).point()[2],
                       (*hit).point()[2] - radius, std::sqrt(
                                ((*hit).point()[0] * (*hit).point()[0]) + ((*hit).point()[1] * (*hit).point()[1]) +
                                ((*hit).point()[2] * (*hit).point()[2])) - radius);
            } else {
                auto point = ray.at(radius);
                nohits_x.push_back(point[0]);
                nohits_y.push_back(point[1]);
                nohits_z.push_back(point[2]);
            }

            Eigen::Array<float, 6, 1> ini;
            ini(Eigen::seq(Eigen::fix<0>, Eigen::fix<2>)) = ray.origin();
            ini(Eigen::seq(Eigen::fix<3>, Eigen::fix<5>)) = ray.direction();
            // Trace curved ray towards de Earth
            for (auto s : method.steps(function, 0.0f, ini, 4.0f * distance)) {
                ray.set_range_max(s.step());

                if (auto hit = earth.trace(ray)) {
                    hits_curved.push_back( Eigen::Array<float,3,1>((*hit).point()[0], (*hit).point()[1], (*hit).point()[2]));
                    hits_nonlinear_x.push_back((*hit).point()[0]);
                    hits_nonlinear_y.push_back((*hit).point()[1]);
                    hits_nonlinear_z.push_back((*hit).point()[2]);
                    printf("Curved hit: X:%f Y:%f Z:%f (Zr:%f)\n", (*hit).point()[0], (*hit).point()[1],
                           (*hit).point()[2], (*hit).point()[2] - radius);
                    break;
                } else {
                    ray = tracer::Ray(s.y()(Eigen::seq(Eigen::fix<0>, Eigen::fix<2>)),
                                      s.y()(Eigen::seq(Eigen::fix<3>, Eigen::fix<5>)));
                }
            }

        }

        /**
        * Figure 7: Shows the footprint of the cone of rays when they intersect with the surface.
        */
        printf("Creating figure...\n");
        // Plot straight hits with the planet (X and Y axis)
        plt_xy.plot(hits_x, hits_y).linestyle("-").color("r").linewidth(1).alpha(alphas[alpha]);
        plt_xy.scatter({center_hit[0]}, {center_hit[1]}).c("r").s(2).alpha(alphas[alpha]);

        // Plot curved hits with the planet (X and Y axis)
        plt_xy.plot(hits_nonlinear_x, hits_nonlinear_y).linestyle("-").color("b").linewidth(1).alpha(alphas[alpha]);
        plt_xy.scatter({center_hit_c[0]}, {center_hit_c[1]}).c("b").s(2).alpha(alphas[alpha]);

        // Plot limits (squared)
        float max_y = *(std::max_element(hits_y.begin(), hits_y.end()));
        float max_x = *(std::max_element(hits_x.begin(), hits_x.end()));
        float min_y = *(std::min_element(hits_y.begin(), hits_y.end()));
        float min_x = *(std::min_element(hits_x.begin(), hits_x.end()));
        printf("Eje X: Desde '%f' hasta '%f'\n", min_x, max_x);
        printf("Eje Y: Desde '%f' hasta '%f'\n", min_y, max_y);

        g_max_x = std::max(max_x, g_max_x);
        g_min_x = std::min(min_x, g_min_x);
        g_max_y = std::max(max_y, g_max_y);
        g_min_y = std::min(min_y, g_min_y);

        alpha += 1;
    }

    float size = std::max(std::max(g_max_x, g_max_y), -std::min(g_min_x, g_min_y));
    std::array<float, 4> limits = {-size, size, -size, size};
    //Plot axis
    plt_xy.plot({-size, size}, {0, 0}).linestyle("-").color("k").linewidth(0.25);
    plt_xy.plot({0, 0}, {-size, size}).linestyle("-").color("k").linewidth(0.25);

    // Add text and save
    plt_xy.xlabel("X");
    plt_xy.ylabel("Y");
    plt_xy.title("Intersection with the Earth surface");
    plt_xy.axis(limits).figsize({477,500}).linewidth(1).savefig("Fig7_3d_cone_graph.svg");

    /** CSV Data save **/
    std::ofstream csv_save;
    csv_save.open ("3d_cone_curved.csv");
    csv_save << "x,y,z\n";
    for (auto hit : hits_curved) {
        csv_save << hit[0] << "," << hit[1] << "," << hit[2] << ",\n";
    }
    csv_save.close();

    csv_save.open ("3d_cone_straight.csv");
    csv_save << "x,y,z\n";
    for (auto hit : hits_straight) {
        csv_save << hit[0] << "," << hit[1] << "," << hit[2] << ",\n";
    }
    csv_save.close();

}