#pragma once
#include <array>

template<typename IOR, typename DIOR>
class Fermat {
    IOR ior;
    DIOR dior;
public:
    Fermat(const IOR& ior, const DIOR& dior) : ior(ior), dior(dior) {}
    
    template<typename Vec6> //First three coordinates x', rest x'' = y'
    Vec6 operator()(const Vec6& v) const {
        Vec6 s;
        for (int i = 0; i<3; ++i) s[i] = v[i+3];
        
        auto n = ior(v[0],v[1],v[2]);
        auto dndr = dior(v[0],v[1],v[2]);
        dndl = dndr[0]*v[4] + dndr[1]*v[5] + dndr[2]*v[6];
        for (int i = 0; i<3; ++i)
            s[i] = (dndr[i] - dndl*v[i+3])/n; 
        return s;
    }
};

template<typename IOR, typename DIOR>
Fermat<IOR,DIOR> fermat(const IOR& ior, const DIOR& dior) {
    return Fermat<IOR,DIOR>(ior,dior);
}

template<typename IOR>
auto fermat(const IOR& ior, float eps = 1.e-5) {
    return fermat(ior, [ior,eps] (float x, float y, float z) 
                { float n = ior(x,y,z); 
                  return std::array<float,3>{(ior(x+eps,y,z) - n)/eps, (ior(x,y+eps,z) - n)/eps, (ior(x,y,z+eps) - n)/eps};
                });
}