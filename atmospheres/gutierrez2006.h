#pragma once

#include "atmosphere.h"

class InversionLayerBase : public pattern::SelfRegisteringReflectableBase {  
public:
    virtual float temperature(float h) const = 0;
    virtual float dtemperature(float h) const = 0;
    static const char* type_name() { return "inversion-layer"; }
    virtual ~InversionLayerBase() = default;
};

class InversionLayer : public pattern::Pimpl<InversionLayerBase> {
public:
    using pattern::Pimpl<InversionLayerBase>::Pimpl;
    using pattern::Pimpl<InversionLayerBase>::operator=;
    
    float temperature(float h) const override {
        return impl()->temperature(h);
    }
    float dtemperature(float h) const override {
        return impl()->dtemperature(h);
    }
};

class InversionLayerFermi : public pattern::SelfRegisteringReflectable<InversionLayerFermi,InversionLayerBase> {
    float height, variation, amplitude;
public:
    InversionLayerFermi(float height = 0, float variation = 30, float amplitude = 20) :
        height(height), variation(variation), amplitude(amplitude) {}
        
    static const char* type_name() { return "fermi"; }
    auto reflect() { return std::tie(height, variation, amplitude); }
    auto reflect_names() const { return std::tuple("height","variation","amplitude"); }
    
    float temperature(float h) const { 
        return -variation+variation/(1.0f+std::exp((h-height)/amplitude));
    }
    
    float dtemperature(float h) const {
        return -(variation*std::exp((h-height)/amplitude))/(amplitude*(1.0f+std::exp((h-height)/amplitude))*(1.0f+std::exp((h-height)/amplitude)));
    }
};

class InversionLayerGaussian : public pattern::SelfRegisteringReflectable<InversionLayerGaussian,InversionLayerBase> {
    float height, variation, amplitude;
public:
    InversionLayerGaussian(float height = 0, float variation = 30, float amplitude = 20) :
        height(height), variation(variation), amplitude(amplitude) {}
        
    static const char* type_name() { return "gaussian"; }
    auto reflect() { return std::tie(height, variation, amplitude); }
    auto reflect_names() const { return std::tuple("height","variation","amplitude"); }
    
    float temperature(float h) const {
        return variation*std::exp(-(h-height)*(h-height)/(200.0*amplitude));
    }
    
    float dtemperature(float h) const {
        return -variation*((h-height)/(100.0*amplitude))*std::exp(-(h-height)*(h-height)/(200.0*amplitude));
    }
};

class Gutierrez2006 : public pattern::SelfRegisteringReflectable<Gutierrez2006,AtmosphereBase>  {
    float planet_radius;
    float wavelength;
    std::list<InversionLayer> inversion_layers;

private:
    static constexpr float lerp(float a, float b, float t) noexcept { return a + t*(b-a);} 
    static constexpr float unlerp(float a, float b, float l) noexcept { return (l-a)/(b-a);} 
    //Data from: https://en.wikipedia.org/wiki/U.S._Standard_Atmosphere
    float temperature_base(float h) const noexcept {
        if (h<0) return 288.15;
        else if (h<11000) return lerp(288.15,216.65,unlerp(0,11000,h));
        else if (h<20000) return lerp(216.65,216.65,unlerp(11000,20000,h));
        else if (h<32000) return lerp(216.65,228.65,unlerp(20000,32000,h));
        else if (h<47000) return lerp(228.65,270.65,unlerp(32000,47000,h));
        else if (h<51000) return lerp(270.65,270.65,unlerp(47000,51000,h));
        else if (h<71000) return lerp(270.65,214.65,unlerp(51000,71000,h));
        else return 214.65;
    }
    float dtemperature_base(float h) const noexcept {
        if (h<0) return 0;
        else if (h<11000) return (216.65-288.15)/11000;
        else if (h<20000) return 0;
        else if (h<32000) return (228.65-216.65)/(32000-20000);
        else if (h<47000) return (270.65-228.65)/(47000-32000);
        else if (h<51000) return 0;
        else if (h<71000) return (214.65-270.65)/(71000-51000);
        else return 0;
    }
    float temperature(float h) const noexcept {
        float t = temperature_base(h);
        for (const InversionLayer& layer : inversion_layers) t+=layer.temperature(h);
        return t;
    }
    float dtemperature(float h) const noexcept {
        float t = dtemperature_base(h);
        for (const InversionLayer& layer : inversion_layers) t+=layer.dtemperature(h);
        return t;
    }
    
    float pressure(float h) const noexcept {
        if (h<0) return 101325;
        else if (h<11000) return lerp(101325,22632.1,unlerp(0,11000,h));
        else if (h<20000) return lerp(22632.1,5474.89,unlerp(11000,20000,h));
        else if (h<32000) return lerp(5474.89,868.019,unlerp(20000,32000,h));
        else if (h<47000) return lerp(868.019,110.906,unlerp(32000,47000,h));
        else if (h<51000) return lerp(110.906,66.9389,unlerp(47000,51000,h));
        else if (h<71000) return lerp(66.9389,3.95642,unlerp(51000,71000,h));
        else return 3.95642;
    }
    float dpressure(float h) const noexcept {
        if (h<0) return 0;
        else if (h<11000) return (22632.1-101325)/11000;
        else if (h<20000) return (5474.89-22632.1)/(20000-11000);
        else if (h<32000) return (868.019-5474.89)/(32000-20000);
        else if (h<47000) return (110.906-868.019)/(47000-32000);
        else if (h<51000) return (66.9389-110.906)/(51000-47000);
        else if (h<71000) return (3.95642-66.9389)/(71000-51000);
        else return 0;
    }    
    
    static constexpr float m = 28.96e-3; 
    static constexpr float r = 8.314510;
    float density(float h) const noexcept {
        return m*pressure(h)/(r*temperature(h));
    }
    float ddensity(float h) const noexcept {
 //       std::cerr<<"h:"<<h<<"  pressure:"<<pressure(h)<<" dpressure="<<dpressure(h)<<" temperature="<<temperature(h)<<" dtemperature="<<dtemperature(h)<<std::endl;
 //       std::cerr<<" ddensity="<<((m/r)*(dpressure(h)*temperature(h) - pressure(h)*dtemperature(h))/(temperature(h)*temperature(h)))<<" Numerical ddensity="<<((density(h+1.e-3)-density(h))/1.e-3)<<std::endl;
 //       std::cin.get();
        return (m/r)*(dpressure(h)*temperature(h) - pressure(h)*dtemperature(h))/(temperature(h)*temperature(h));
    }
    
    float ior_base() const noexcept {
        //Cauchy's dispersion formula from the paper
        //const float a = 28.79e-5; const float b = 5.67e-5;
        //return a*(1.0f+b/(wavelength*wavelength)) + 1.0f;
        //Corrected Cauchy formula with data from https://hypertextbook.com/facts/2005/MayaBarsky.shtml
        const float b = 0.00004989*324*1.e-14/77;
        const float a = 1.00032408 - b*0.25e14; 
        return a+b/(wavelength*wavelength);
    }
    float ior_height(float h) const noexcept {
        //Gladstone-Dale
        return density(h)*(ior_base()-1.0f)+1.0f;
    } 
    float dior_height(float h) const noexcept {
        //Gladstone-Dale
        return ddensity(h)*(ior_base()-1.0f);       
    }     
    

public:
    Gutierrez2006(const float planet_radius = 6370949) : planet_radius(planet_radius),wavelength(500.e-9) {}

    float ior (float x, float y, float z) const override { // Index of Refraction
        float h = std::sqrt(x*x + y*y + z*z) - planet_radius;
        return ior_height(h); // 1.0+1.e-6 for refractivity->IOR
    };

    std::array<float, 3> dior (float x, float y, float z) const override { // IOR gradient
        float h = std::sqrt(x*x + y*y + z*z) - planet_radius;
        float k = dior_height(h)/std::sqrt(x*x+y*y+z*z); // 1.e-6 factor for drefractivity->dIOR
   //     std::cerr<<x<<", "<<y<<", "<<z<<" -> "<<k*x<<", "<<k*y<<", "<<k*z<<" vs "<<((ior(x+1.e-1,y,z)-ior(x,y,z))/1.e-1)<<", "<<((ior(x,y+1.e-1,z)-ior(x,y,z))/1.e-1)<<", "<<((ior(x,y,z+1.e-1)-ior(x,y,z))/1.e-1)<<std::endl;

        return std::array<float,3>{k*x,k*y,k*z};
    };
    
    static const char* type_name() { return "gutierrez2006"; }    
    auto reflect_names() const { return std::tuple("wavelength"); }
    auto reflect() { return std::tie(wavelength,inversion_layers); }


};



