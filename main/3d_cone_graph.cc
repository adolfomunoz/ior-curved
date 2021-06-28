#include <svg-cpp-plot/svg-cpp-plot.h>
#include <ivp/ivp.h>
#include <mj2/tracer/primitives/sphere.h>
#include <mj2/tracer/primitives/plane.h>
#include "../eq/fermat.h"
#include <cmath>
#include <string>
/**
 * Plot [Figure 7] that shows in overhead view (XY axis) the hits of the Cherenkov cone using the three dimensions.
 */
int main(int argc, char** argv) {
    /** Parameters */
    // Plot
    //float width = 200;
    //float height = 50;

    // Rays
    double angle = 0*(M_PI/180.0); // Zenith angle to radians
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
    double omega_ch = 1*M_PI/180.0;
    for (int i = 0; i<(argc-1); ++i) { // Custom angle as argument
        if (std::string(argv[i]) == "-omega_ch") omega_ch = atof(argv[++i])*M_PI/180.0;
    }
    printf("omega_ch: %f degrees\n", omega_ch*(180/M_PI));

    // Station basis
    Eigen::Vector3f o = Eigen::Vector3f(origin[0], origin[1], origin[2]);
    Eigen::Vector3f i = Eigen::Vector3f(-std::sin(angle), 0, -std::cos(angle));
    i = i.normalized();
    Eigen::Vector3f j = Eigen::Vector3f(1, 1, (i[0] + i[1])/-i[2]);
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
    m <<    i[0], j[0], k[0],
            i[1], j[1], k[1],
            i[2], j[2], k[2];

    // Center ray towards Earth
    Eigen::Vector3f dir_c = m * Eigen::Vector3f(1, 0, 0);
    printf("Center ray direction: %f, %f, %f\n", dir_c[0], dir_c[1], dir_c[2]);
    tracer::Ray center_ray(origin, dir_c);
    if (auto hit = earth.trace(center_ray)) {
        printf("Center Hit: X:%f Y:%f Z:%f\n", (*hit).point()[0], (*hit).point()[1], (*hit).point()[2]);
        printf("Center Hit (z=z-radius): X:%f Y:%f Z:%f\n", (*hit).point()[0], (*hit).point()[1], (*hit).point()[2]-radius);
    }

    std::list<float> zeniths, angles, omegas, omegas_c, rhos, rhos_c, rhos_diff;
    std::list<float> hits_x, hits_y, hits_z, hits_nonlinear_x, hits_nonlinear_y, hits_nonlinear_z;
    std::list<float> xc,yc,zc;
    for (float c = 0; c < 2*M_PI; c += M_PI_4/16) {
        // Cone calculation
        Eigen::Vector3f v = Eigen::Vector3f(std::cos(omega_ch),std::sin(omega_ch)*std::cos(c),std::sin(omega_ch)*std::sin(c));
        //printf("v: %f, %f, %f\n", v[0], v[1], v[2]);
        Eigen::Vector3f dir = m * v;
        xc.push_back(dir[0]);
        yc.push_back(dir[1]);
        zc.push_back(dir[2]);
        //printf("dir: %f, %f, %f\n", dir[0], dir[1], dir[2]);

        // Trace ray towards the Earth
        tracer::Ray ray(origin, dir);
        if (auto hit = earth.trace(ray)) {
            hits_x.push_back((*hit).point()[0]);
            hits_y.push_back((*hit).point()[1]);
            hits_z.push_back((*hit).point()[2]);
            printf("Hit: X:%f Y:%f Z:%f (Zr:%f) (%f)\n", (*hit).point()[0], (*hit).point()[1], (*hit).point()[2], (*hit).point()[2]-radius, std::sqrt(((*hit).point()[0]*(*hit).point()[0])+((*hit).point()[1]*(*hit).point()[1])+((*hit).point()[2]*(*hit).point()[2]))-radius);
        }

        Eigen::Array<float, 6, 1> ini;
        ini(Eigen::seq(Eigen::fix<0>, Eigen::fix<2>)) = ray.origin();
        ini(Eigen::seq(Eigen::fix<3>, Eigen::fix<5>)) = ray.direction();
        // Trace curved ray towards de Earth
        for (auto s : method.steps(function, 0.0f, ini, 4.0f * distance)) {//1.0f*float(atmosphere_height))) {
            ray.set_range_max(s.step());

            if (auto hit = earth.trace(ray)) {
                hits_nonlinear_x.push_back((*hit).point()[0]);
                hits_nonlinear_y.push_back((*hit).point()[1]);
                printf("Curved hit: X:%f Y:%f Z:%f (Zr:%f)\n", (*hit).point()[0], (*hit).point()[1], (*hit).point()[2], (*hit).point()[2]-radius);
                break;
            } else {
                //printf("No hit: X:%f Y:%f Z:%f (Zr:%f) (%f)\n", s.y()[0], s.y()[1], s.y()[2], s.y()[2]-radius, std::sqrt((s.y()[0]*s.y()[0])+(s.y()[1]*s.y()[1])+(s.y()[2]*s.y()[2]))-radius);
            }
        }
    }

    /**
    * Figure 7: Shows the
    */
    printf("Creating figure...\n");
    svg_cpp_plot::SVGPlot plt_xy;
    // Use relative to the maximum Z
    /*float max_z = *(std::max_element(hits_z.begin(), hits_z.end()));
    std::list<float> hits_zz;
    for (auto &zz : hits_z) {
        hits_zz.push_back(-(max_z-zz));
    }*/
    // Plot straight hits with the planet (X and Y axis)
    plt_xy.plot(hits_x, hits_y).linestyle("-").color( "r").linewidth(1);
    plt_xy.scatter(hits_x,hits_y).c("r").s(2).alpha(0.5);
    // Plot curved hits with the planet (X and Y axis)
    plt_xy.plot(hits_nonlinear_x, hits_nonlinear_y).linestyle("-").color( "b").linewidth(1);
    plt_xy.scatter(hits_nonlinear_x,hits_nonlinear_y).c("b").s(2).alpha(0.5);
    // Plot limits (squared)
    float max_y = *(std::max_element(hits_y.begin(),hits_y.end()));
    float max_x = *(std::max_element(hits_x.begin(), hits_x.end()));
    float min_y = *(std::min_element(hits_y.begin(),hits_y.end()));
    float min_x = *(std::min_element(hits_x.begin(), hits_x.end()));
    printf("Eje X: Desde '%f' hasta '%f'\n",min_x,max_x);
    printf("Eje Y: Desde '%f' hasta '%f'\n",min_y,max_y);
    float size = std::max(std::max(max_x,max_y),-std::min(min_x,min_y));
    std::array<float,4> limits = {-size,size,-size,size};
    //Plot axis
    plt_xy.plot({-size,size}, {0,0}).linestyle("-").color( "k").linewidth(0.25);
    plt_xy.plot({0,0}, {-size,size}).linestyle("-").color( "k").linewidth(0.25);
    // Add text and save
    //plt_zenith.xlabel("X");
    //plt_zenith.ylabel("Y");
    //plt_zenith.title("Circle");
    //plt_zenith.linewidth(1).savefig("3d_cone_graph.svg");
    //plt_zenith.xticks({}).yticks({}).figsize({500,500}).linewidth(0).savefig("3d_cone_graph.svg");
    plt_xy.axis(limits).xticks({}).yticks({}).figsize({500,500}).linewidth(1).savefig("Fig7_3d_cone_graph.svg");

}