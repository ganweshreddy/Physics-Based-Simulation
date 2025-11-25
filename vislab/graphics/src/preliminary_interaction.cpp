#include <vislab/graphics/preliminary_interaction.hpp>

#include <vislab/graphics/ray.hpp>
#include <vislab/graphics/shape.hpp>

namespace vislab
{
    bool PreliminaryIntersection::isValid() const { return t != std::numeric_limits<double>::infinity(); }

    SurfaceInteraction PreliminaryIntersection::computeSurfaceInteraction(const Ray3d& ray, EHitComputeFlag flags)
    {
        const Shape* target   = shape;
        SurfaceInteraction si = target->computeSurfaceInteraction(ray, *this, flags);
        bool active           = si.isValid();

        if (!active)
            si.t = std::numeric_limits<double>::infinity();
        si.prim_index = prim_index;

        // Set shape pointer if not already set by compute_surface_interaction()
        if (si.shape == nullptr)
            si.shape = shape;

        if (hasFlag(flags, EHitComputeFlag::ShadingFrame))
            si.initializeShadingFrame();

        // Incident direction in local coordinates
        si.wi = active ? si.toLocal(-ray.direction) : -ray.direction;

        si.duv_dx = Eigen::Vector2d::Zero();
        si.duv_dy = Eigen::Vector2d::Zero();
        return si;
    }
}
