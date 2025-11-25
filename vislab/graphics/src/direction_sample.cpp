#include <vislab/graphics/direction_sample.hpp>

#include <vislab/graphics/ray.hpp>
#include <vislab/graphics/surface_interaction.hpp>

namespace vislab
{

    DirectionSample::DirectionSample()
        : PositionSample()
        , direction(0, 0, 0)
        , dist(0)
        , light(nullptr)
    {
    }

    DirectionSample::DirectionSample(const Scene* scene, const SurfaceInteraction& si, const Interaction& it)
        : PositionSample(si)
    {
        direction = si.position - it.position;
        direction.normalize();
        if (!si.isValid())
            direction = -si.wi;
        light = si.light(scene);
    }

    DirectionSample::DirectionSample(const Eigen::Vector3d& p, const Eigen::Vector3d& n, const Eigen::Vector2d& uv, const double& pdf, const bool& delta, const Eigen::Vector3d& d, const double& dist, const Light* light)
        : PositionSample(p, n, uv, pdf, delta)
        , direction(d)
        , dist(dist)
        , light(light)
    {
    }

    DirectionSample::DirectionSample(const PositionSample& base)
        : PositionSample(base)
        , direction(0, 0, 0)
        , dist(0)
        , light(nullptr)
    {
    }
}
