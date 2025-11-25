#include <vislab/graphics/position_sample.hpp>

#include <vislab/graphics/surface_interaction.hpp>

namespace vislab
{
    PositionSample::PositionSample()
        : position(0, 0, 0)
        , normal(0, 0, 0)
        , uv(0, 0)
        , pdf(0)
        , delta(false)
    {
    }

    PositionSample::PositionSample(const SurfaceInteraction& si)
        : position(si.position)
        , normal(si.sh_frame.n)
        , uv(si.uv)
        , pdf(0)
        , delta(false)
    {
    }

    PositionSample::PositionSample(const Eigen::Vector3d& p, const Eigen::Vector3d& n, const Eigen::Vector2d& uv, const double& pdf, const bool& delta)
        : position(p)
        , normal(n)
        , uv(uv)
        , pdf(pdf)
        , delta(delta)
    {
    }
}
