#include <vislab/graphics/area_light.hpp>

#include <vislab/graphics/interaction.hpp>
#include <vislab/graphics/ray.hpp>
#include <vislab/graphics/shape.hpp>
#include <vislab/graphics/surface_interaction.hpp>
#include <vislab/graphics/warp.hpp>

namespace vislab
{
    AreaLight::AreaLight()
        : radiance(Spectrum::Zero())
    {
    }

    std::pair<Ray3d, Spectrum> AreaLight::sampleRay(const Eigen::Vector2d& sample, const Eigen::Vector2d& sample2) const
    {
        PositionSample ps = shape->samplePosition(sample);
        SurfaceInteraction si(ps);

        Eigen::Vector3d local = Warp::square_to_cosine_hemisphere(sample);
        return std::make_pair(
            Ray3d(si.position, si.toWorld(local)),
            radiance * EIGEN_PI / ps.pdf);
    }

    std::pair<DirectionSample, Spectrum> AreaLight::sampleDirection(const Interaction& it, const Eigen::Vector2d& sample) const
    {
        DirectionSample ds = shape->sampleDirection(it, sample);
        if (ds.direction.dot(ds.normal) >= 0. || ds.pdf == 0.)
            return { ds, Spectrum::Zero() };

        SurfaceInteraction si(ds);
        Spectrum spec = radiance / ds.pdf;
        ds.light      = this;
        return { ds, spec };
    }

    double AreaLight::pdfDirection(const Interaction& it, const DirectionSample& ds) const
    {
        double dp = ds.direction.dot(ds.normal);
        if (dp < 0.)
            return 0.;
        return shape->pdfDirection(it, ds);
    }

    Spectrum AreaLight::evaluate(const SurfaceInteraction& si) const
    {
        return Frame::cos_theta(si.wi) > 0.f ? radiance : Spectrum::Zero();
    }

    uint32_t AreaLight::flags() const
    {
        return (uint32_t)ELightFlag::Surface;
    }
}
