#include <vislab/graphics/light.hpp>

namespace vislab
{
    bool Light::isDeltaPosition() const { return hasFlag(ELightFlag::DeltaPosition); }
    bool Light::isDeltaDirection() const { return hasFlag(ELightFlag::DeltaDirection); }
    bool Light::isArea() const { return hasFlag(ELightFlag::Surface); }
    bool Light::isInfinite() const { return hasFlag(ELightFlag::Infinite); }
    bool Light::isDeltaLight() const { return isDeltaPosition() || isDeltaDirection(); }

    bool Light::isEnvironment() const
    {
        return hasFlag(ELightFlag::Infinite) && !hasFlag(ELightFlag::Delta);
    }

    bool Light::hasFlag(const ELightFlag& flag) const
    {
        uint32_t f = (uint32_t)flag;
        return (flags() & f) != 0;
    }
}
