#include <vislab/graphics/point_light.hpp>

#include <vislab/graphics/interaction.hpp>
#include <vislab/graphics/ray.hpp>
#include <vislab/graphics/surface_interaction.hpp>
#include <vislab/graphics/warp.hpp>

namespace vislab
{
    PointLight::PointLight()
        : intensity(Eigen::Vector3d::Zero())
    {
    }

    std::pair<Ray3d, Spectrum> PointLight::sampleRay(const Eigen::Vector2d& sample, const Eigen::Vector2d& sample2) const
    {
        Eigen::Vector3d origin    = transform.getMatrix().translation();
        Eigen::Vector3d direction = Warp::square_to_uniform_sphere(sample);
        Spectrum power            = intensity * 4. * EIGEN_PI;
        return { Ray3d(origin, direction), power };
    }

    std::pair<DirectionSample, Spectrum> PointLight::sampleDirection(const Interaction& it, const Eigen::Vector2d& sample) const
    {
        DirectionSample ds;
        ds.position  = transform.getMatrix().translation();
        ds.normal    = Eigen::Vector3d::Zero();
        ds.uv        = Eigen::Vector2d::Zero();
        ds.pdf       = 1.f;
        ds.delta     = true;
        ds.light     = this;
        ds.direction = ds.position - it.position;

        double dist2    = ds.direction.squaredNorm(),
               inv_dist = 1. / std::sqrt(dist2);

        ds.dist = std::sqrt(dist2);
        ds.direction *= inv_dist;

        Spectrum spec = intensity * inv_dist * inv_dist;
        return { ds, spec };
    }

    double PointLight::pdfDirection(const Interaction& it, const DirectionSample& ds) const
    {
        return 0.f;
    }

    Spectrum PointLight::evaluate(const SurfaceInteraction& si) const
    {
        return Spectrum::Zero();
    }

    uint32_t PointLight::flags() const
    {
        return (uint32_t)ELightFlag::DeltaPosition;
    }
}
