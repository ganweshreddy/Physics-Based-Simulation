#include <vislab/graphics/surface_interaction.hpp>

#include <vislab/graphics/bsdf.hpp>
#include <vislab/graphics/position_sample.hpp>
#include <vislab/graphics/ray_differential.hpp>
#include <vislab/graphics/scene.hpp>
#include <vislab/graphics/shape.hpp>

namespace vislab
{
    SurfaceInteraction::SurfaceInteraction()
        : Interaction()
        , wi(0, 0, 0)
        , normal(0, 0, 0)
        , shape(NULL)
        , prim_index(0)
        , uv(0, 0)
        , sh_frame()
        , dp_du(0, 0, 0)
        , dp_dv(0, 0, 0)
        , duv_dx(0, 0)
        , duv_dy(0, 0)
    {
    }

    SurfaceInteraction::SurfaceInteraction(const PositionSample& ps)
        : Interaction(0, ps.position)
        , uv(ps.uv)
        , normal(ps.normal)
        , sh_frame(ps.normal)
        , prim_index(0)
        , shape(NULL)
        , dp_du(0, 0, 0)
        , dp_dv(0, 0, 0)
        , duv_dx(0, 0)
        , duv_dy(0, 0)
    {
    }

    Eigen::Vector3d SurfaceInteraction::toWorld(const Eigen::Vector3d& v) const
    {
        return sh_frame.toWorld(v);
    }

    Eigen::Vector3d SurfaceInteraction::toLocal(const Eigen::Vector3d& v) const
    {
        return sh_frame.toLocal(v);
    }

    const Light* SurfaceInteraction::light(const Scene* scene) const
    {
        if (isValid())
            return shape->light.get();
        else
            return scene->environment.get();
    }

    void SurfaceInteraction::initializeShadingFrame()
    {
        sh_frame.s = (-sh_frame.n * sh_frame.n.dot(dp_du) + dp_du).normalized();
        sh_frame.t = sh_frame.n.cross(sh_frame.s);
    }

    const BSDF* SurfaceInteraction::bsdf(const RayDifferential3d& ray)
    {
        const BSDF* bsdf = shape->bsdf.get();

        if (!hasUVPartials() && bsdf->needsDifferentials())
            computeUVPartials(ray);

        return bsdf;
    }

    const BSDF* SurfaceInteraction::bsdf() const
    {
        return shape->bsdf.get();
    }

    void SurfaceInteraction::computeUVPartials(const RayDifferential3d& ray)
    {
        if (!ray.hasDifferentials)
            return;

        // Compute interaction with the two offset rays
        double d   = normal.dot(position),
               t_x = (d - normal.dot(ray.origin_x)) / normal.dot(ray.direction_x),
               t_y = (d - normal.dot(ray.origin_y)) / normal.dot(ray.direction_y);

        // Corresponding positions near the surface
        Eigen::Vector3d dp_dx = ray.direction_x * t_x + ray.origin_x - position,
                        dp_dy = ray.direction_y * t_y + ray.origin_y - position;

        // Solve a least squares problem to turn this into UV coordinates
        double a00     = dp_du.dot(dp_du),
               a01     = dp_du.dot(dp_dv),
               a11     = dp_dv.dot(dp_dv),
               inv_det = 1.0 / (a00 * a11 - a01 * a01);

        double b0x = dp_du.dot(dp_dx),
               b1x = dp_dv.dot(dp_dx),
               b0y = dp_du.dot(dp_dy),
               b1y = dp_dv.dot(dp_dy);

        // Set the UV partials to zero if dpdu and/or dpdv == 0
        inv_det = inv_det == inv_det ? inv_det : 0.f;

        duv_dx = Eigen::Vector2d(a11 * b0x - a01 * b1x,
                                 a00 * b1x - a01 * b0x) *
                 inv_det;

        duv_dy = Eigen::Vector2d(a11 * b0y - a01 * b1y,
                                 a00 * b1y - a01 * b0y) *
                 inv_det;
    }

    bool SurfaceInteraction::hasUVPartials() const
    {
        return duv_dx.any() || duv_dy.any();
    }

    Ray3d SurfaceInteraction::spawnRay(const Eigen::Vector3d& direction) const
    {
        Ray3d ray;
        ray.origin    = offsetPosition(direction);
        ray.direction = direction;
        ray.tMax      = std::numeric_limits<double>::max();
        return ray;
    }

    Eigen::Vector3d SurfaceInteraction::offsetPosition(const Eigen::Vector3d& direction) const
    {
        double mag   = (1. + position.cwiseAbs().maxCoeff()) * rayEpsilon();
        double NdotD = normal.dot(direction);
        if (NdotD < 0)
            mag *= -1;
        return mag * normal + position;
    }
}
