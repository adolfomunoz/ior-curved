#pragma once

#include <patterns/patterns.h>
#include <array>

class AtmosphereBase : public pattern::SelfRegisteringReflectableBase {
public:
    static const char* type_name() { return "atmosphere"; }
    virtual float ior (float x, float y, float z) const = 0;
    virtual std::array<float, 3> dior (float x, float y, float z) const = 0;
    virtual ~AtmosphereBase() = default;

};

class Atmosphere : public pattern::Pimpl<AtmosphereBase> {
public:
    using pattern::Pimpl<AtmosphereBase>::Pimpl;
    using pattern::Pimpl<AtmosphereBase>::operator=;
    
    float ior (float x, float y, float z) const override {
        return impl()->ior(x,y,z);
    }
    std::array<float, 3>  dior (float x, float y, float z) const override {
        return impl()->dior(x,y,z);
    }
};


