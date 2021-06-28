//
// Created by workccu on 28/06/2021.
//

#include <svg-cpp-plot/svg-cpp-plot.h>
#include <ivp/ivp.h>
#include <mj2/tracer/primitives/sphere.h>
#include "../eq/fermat.h"
#include <cmath>
#include <string>
// Tested values: -zenith 88 -omega_ch 1
int main(int argc, char** argv) {
    /** Parameters */
    // Plot
    float width = 200;
    float height = 200;

    // Rays
    double angle = 88*M_PI/180.0; // Zenith angle to radians
    for (int i = 0; i<(argc-1); ++i) { // Custom angle as argument
        if (std::string(argv[i]) == "-zenith") angle = atof(argv[++i])*M_PI/180.0;
    }
    printf("Zenith: %f degrees\n", angle*(180/M_PI));

    // Scene
    const double radius = 6370949; // Earth radius (m)
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
    float distance = (-10 + 12 / std::cos(angle)) * 1000; // x1000 km -> m
    float x = distance * std::sin(angle);
    float y = 0;
    float z = ( distance * std::cos(angle) ) + radius;
    float atmosphere_height = std::sqrt(x*x + z*z) - radius;

    printf("Origin of rays: %f,%f\n", x, z);
    Eigen::Vector3f origin(x,0,z);
    printf("Distance Origin-North: %f m\n", distance);
    printf("Atmosphere height: %f m\n", atmosphere_height);

    auto function = fermat(ior,dior);
    //IVP::Adaptive<IVP::Dopri> method(100,1.e-3); // (Pasos inicialmente, tolerancia)
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

    svg_cpp_plot::SVGPlot pltxy, pltxz;
    std::list<float> zeniths, angles, omegas, omegas_c, rhos, rhos_c, rhos_diff;
    std::list<float> hits_x, hits_y, hits_z, hits_nonlinear_x, hits_nonlinear_y, hits_nonlinear_z, nohits_x, nohits_y, path_x, path_y, path_z, nohits_z;
    std::list<float> xc,yc,zc;
    for (float c = 0; c < 2*M_PI; c += M_PI/4) {
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
        std::cout << "Straight ray origin:\n" << ray.origin() << "\n";
        std::cout << "Straight ray direction:\n" << ray.direction() << "\n";
        if (auto hit = earth.trace(ray)) {
            hits_x.push_back((*hit).point()[0]);
            hits_y.push_back((*hit).point()[1]);
            hits_z.push_back((*hit).point()[2]);
            printf("Hit: X:%f Y:%f Z:%f (Zr:%f) (%f)\n", (*hit).point()[0], (*hit).point()[1], (*hit).point()[2], (*hit).point()[2]-radius, std::sqrt(((*hit).point()[0]*(*hit).point()[0])+((*hit).point()[1]*(*hit).point()[1])+((*hit).point()[2]*(*hit).point()[2]))-radius);
        } else {
            auto point = ray.at(radius);
            nohits_x.push_back(point[0]);
            nohits_y.push_back(point[1]);
            nohits_z.push_back(point[2]);
        }

        Eigen::Array<float, 6, 1> ini;
        ini(Eigen::seq(Eigen::fix<0>, Eigen::fix<2>)) = ray.origin();
        ini(Eigen::seq(Eigen::fix<3>, Eigen::fix<5>)) = ray.direction();
        std::cout << "Curved ray origin:\n" << ray.origin() << "\n";
        // Trace curved ray towards de Earth
        bool debug = false;
        for (auto s : method.steps(function, 0.0f, ini, 4.0f * distance)) {//1.0f*float(atmosphere_height))) {
            ray.set_range_max(s.step());

            if (!debug) {
                //std::cout << "Curved ray direction:\n" << ray.direction() << "\n";
                debug = false;
            }
            if (auto hit = earth.trace(ray)) {
                hits_nonlinear_x.push_back((*hit).point()[0]);
                hits_nonlinear_y.push_back((*hit).point()[1]);
                hits_nonlinear_z.push_back((*hit).point()[2]);
                printf("Curved hit: X:%f Y:%f Z:%f (Zr:%f)\n", (*hit).point()[0], (*hit).point()[1], (*hit).point()[2], (*hit).point()[2]-radius);
                path_x.push_back(hits_nonlinear_x.back());
                path_y.push_back(hits_nonlinear_y.back());
                path_z.push_back(hits_nonlinear_z.back());
                break;
            } else {
                //std::cout << "Debug:\n" <<s.y() << "\n";
                //std::cout << "Curved ray direction:\n" << s.y()[0] << " " << s.y()[1] << " " << s.y()[2]-radius << " | " << s.y()[3] << " " << s.y()[4] << " " << s.y()[5] << "\n";
                path_x.push_back(s.y()[0]);
                path_y.push_back(s.y()[1]);
                path_z.push_back(s.y()[2]);
                //printf("No hit: X:%f Y:%f Z:%f (Zr:%f) (%f)\n", s.y()[0], s.y()[1], s.y()[2], s.y()[2]-radius, std::sqrt((s.y()[0]*s.y()[0])+(s.y()[1]*s.y()[1])+(s.y()[2]*s.y()[2]))-radius);

            }
        }
        // Plot curved ray from origin towards Earth
        //pltxy.plot(path_x,path_y).color("r").linewidth(0.25);
        //pltxz.plot(path_x,path_z).color("r").linewidth(0.25);
        pltxy.scatter(path_x,path_y).c("r").s(0.2); // Using scatters because plot causing unwanted "return to origin" lines. Need to see this TODO
        pltxz.scatter(path_x,path_z).c("r").s(0.2);

    }


    // Plot Earth
    std::list<float> earth_x, earth_y;
    for (float a = -std::asin(x/radius); a<=std::asin(x/radius); a+=M_PI/180) {
        earth_x.push_back(radius*std::sin(a));
        earth_y.push_back(radius*std::cos(a));
    }
    pltxy.plot(earth_x,earth_y).color("b").linewidth(0.1);

    // Get the limits of the plot
    float min_x = *(std::min_element(hits_x.begin(),hits_x.end()));
    float min_y = *(std::min_element(hits_y.begin(),hits_y.end()));
    std::array<float,4> limits{min_x,x,min_y,y};
    if ((limits[1]-limits[0])*width>(limits[3]-limits[2])*height)
        limits[2] = limits[3] - (limits[1]-limits[0])*height/width;
    else
        limits[0] = limits[1] - (limits[3]-limits[2])*width/height;

    float d = (limits[3]-limits[2])/32.0;
    limits[0] -= d; limits[1] += d; limits[2] -=d; limits[3] += d;


    // Plot rays from origin towards Earth as discontinued line
    std::list<float>::const_iterator ie,je;
    for (ie = hits_x.begin(), je = hits_y.begin(); (ie != hits_x.end()) && (je != hits_y.end()); ++ie, ++je)
        pltxy.plot({x,*ie},{y,*je}).color("k").linewidth(0.25).format("--");
    for (ie = nohits_x.begin(), je = nohits_y.begin(); (ie != nohits_x.end()) && (je != nohits_y.end()); ++ie, ++je)
        pltxy.plot({x,*ie},{y,*je}).color("k").linewidth(0.25).format("--");

    pltxy.scatter(hits_x,hits_y).c("k"); // Plot point where ray intersect Earth
    pltxy.scatter(hits_nonlinear_x,hits_nonlinear_y).c("r"); // Plot point where curved ray intersect Earth
    pltxy.scatter({x}, {y}).c("g"); // Plot origin of the rays


    // Create figure and save it.
    printf("Creating figure...");
    pltxy.axis(limits).linewidth(0).xticks({}).yticks({}).figsize({width,height}).savefig("atmosphere_dxy.svg");





    // Plot Earth
    for (float a = -std::asin(x/radius); a<=std::asin(x/radius); a+=M_PI/180) {
        earth_x.push_back(radius*std::sin(a));
        earth_y.push_back(radius*std::cos(a));
    }
    pltxz.plot(earth_x,earth_y).color("b").linewidth(0.1);

    // Get the limits of the plot
    min_x = *(std::min_element(hits_x.begin(),hits_x.end()));
    min_y = *(std::min_element(hits_z.begin(),hits_z.end()));
    limits = {min_x,x,min_y,z};
    if ((limits[1]-limits[0])*width>(limits[3]-limits[2])*height)
        limits[2] = limits[3] - (limits[1]-limits[0])*height/width;
    else
        limits[0] = limits[1] - (limits[3]-limits[2])*width/height;

    d = (limits[3]-limits[2])/32.0;
    limits[0] -= d; limits[1] += d; limits[2] -=d; limits[3] += d;


    // Plot rays from origin towards Earth as discontinued line

    for (ie = hits_x.begin(), je = hits_z.begin(); (ie != hits_x.end()) && (je != hits_z.end()); ++ie, ++je)
        pltxz.plot({x,*ie},{z,*je}).color("k").linewidth(0.25).format("--");
    for (ie = nohits_x.begin(), je = nohits_z.begin(); (ie != nohits_x.end()) && (je != nohits_z.end()); ++ie, ++je)
        pltxz.plot({x,*ie},{z,*je}).color("k").linewidth(0.25).format("--");

    pltxz.scatter(hits_x,hits_z).c("k"); // Plot point where ray intersect Earth
    pltxz.scatter(hits_nonlinear_x,hits_nonlinear_z).c("r"); // Plot point where curved ray intersect Earth
    pltxz.scatter({x}, {z}).c("g"); // Plot origin of the rays


    // Create figure and save it.
    printf("Creating figure...");
    pltxz.axis(limits).linewidth(0).xticks({}).yticks({}).figsize({width,height}).savefig("atmosphere_dxz.svg");

}

