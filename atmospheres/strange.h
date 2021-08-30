class Strange {
    float planet_radius;

private:
    // taking data from https://www.researchgate.net/publication/305788876_The_Tropical_Air-Sea_Propagation_Study_TAPS
    // this is refractivity over sea, and is an extreme case. Im parametrizing the data in figure 3 with 4 lines passing on points (Altitude,M):
    //  (0,375) -> (275,412.5)
    //  (275,412.5) -> (350,416)
    //  (350,416) -> (425,384)
    //  (425,384) -> (600,412.5)
    //  M=N+0.157*h
    //   y-y1=((y2-y1)/(x2-x1))(x-x1)
    // (otra para mirar si seguimos por aca:Analysis of the Target Detection Performance of Air-to-Air Airborne Radar Using Long-Range3Propagation Simulation in Abnormal Atmospheric Conditions)

    // h in meters
    float refractivity(float h) {
        float a, M, N = 0;
        if (h<275.0) {
            // p1=(0,375) p2=(275,412.5)
            a=(412.5-375.0)/(275.0-0.0);
            M= a*h - a*0.0 + 375.0;
            N= M - 0.157*h;
            return N;
        }
        if (h<350.0) {
            // (275,412.5) -> (350,416)
            a=(416.0-412.5)/(350.0-275.0);
            M= a*h - a*275.0 + 412.5;
            N= M - 0.157*h;
            return N;
        }
        if(h<425.0) {
            //    (350,416.0) -> (425,384)
            a=(384.0-416.0)/(425.0-350.0);
            M= a*h - a*350.0 + 416.0;
            N= M - 0.157*h;
            return N ;
        }
        if(h<600.0) {
            //  (425,384) -> (600,412.5)
            a=(412.5-384.0)/(600.0-425.0);
            M= a*h - a*425.0 + 384.0;
            N= M - 0.157*h;
            return N ;
        }
        //  i will use N0 as 375 to make it coincede more or less with the  measurement, and keep k as 0.1218 to change as little as possible, but
        // if N=375exp(-k h), at 600m this has to mate with 412,5-0.157*600, -> 318,3. but the standard exponential model gives 348,5724, so ill displace h, 318.3=375.exp(-k (600+d))  -ln(318.3/375)/k -600 = d
        N=375.0*std::exp(-0.0001218 *(h+745.908798706));
        return N;
    }

    // h in meters
    float refractivityDerivative(float h) {
        float a, N;
        if(h<275.0) {
            // p1=(0,375) p2=(275,412.5)
            a=(412.5-375.0)/(275.0-0.0);
            N= a - 0.157;
            return N;
        }
        if(h<350.0) {
            // (275,416) -> (350,412.5)
            //    N= - 0.157
            a=(416.0-412.5)/(350.0-275.0);
            N= a - 0.157;
            return N;
        }
        if(h<425.0) {
            //    (350,416) -> (425,384)
            a=(384.0-416.)/(425.0-350.0);
            N= a - 0.157;
            return N;
        }
        if(h<600.0) {
            //  (425,384) -> (600,412.5)
            a=(412.5-384.0)/(600.0-425.0);
            N= a - 0.157;
            return N;
        }
        // if N=375exp(-k h), at 600m this has to mate with 412,5-0.157*600, -> 318,3. but the standard exponential model gives 348,5724, so ill displace h, 318.3=375.exp(-k (600+d))  -ln(318.3/375)/k -600 = d
        N=-0.0001218*375.0*std::exp(-0.0001218 *(h+745.908798706));
        return N;
    }

public:
    Strange(const float planet_radius) : planet_radius(planet_radius) {}

    float ior (float x, float y, float z) { // Index of Refraction
        float h = std::sqrt(x*x + y*y + z*z) - planet_radius;
        return 1.0+1.e-6*refractivity(h); // 1.0+1.e-6 for refractivity->IOR
    };

    std::array<float, 3> dior (float x, float y, float z) { // IOR gradient
        float h = std::sqrt(x*x + y*y + z*z) - planet_radius;
        float k = (1.e-6*refractivityDerivative(h))/std::sqrt(x*x+y*y+z*z); // 1.e-6 factor for drefractivity->dIOR
        return std::array<float,3>{k*x,k*y,k*z};
    };

};



