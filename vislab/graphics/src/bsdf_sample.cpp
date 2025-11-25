#include <vislab/graphics/bsdf_sample.hpp>

namespace vislab
{
    bool BSDFSample::hasFlag(const EBSDFFlag& flag) const
    {
        uint32_t f = (uint32_t)flag;
        return (sampledType & f) != 0;
    }

    BSDFSample::BSDFSample(const Eigen::Vector3d& wo)
        : wo(wo)
        , pdf(0.)
        , eta(1.)
        , sampledType(0)
    {
    }
}
