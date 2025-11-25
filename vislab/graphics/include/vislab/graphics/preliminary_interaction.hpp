#pragma once

#include "hit_compute_flag.hpp"
#include "ray_fwd.hpp"
#include "surface_interaction.hpp"

namespace vislab
{
    class Shape;

    /**
     * @brief Stores preliminary information related to a ray intersection.
     */
    struct PreliminaryIntersection
    {
        /**
         * @brief Distance traveled along the ray.
         */
        double t = std::numeric_limits<double>::infinity();

        /**
         * @brief 2D coordinates on the primitive surface parameterization.
         */
        Eigen::Vector2d prim_uv;

        /**
         * @brief Primitive index, e.g. the triangle ID (if applicable).
         */
        uint32_t prim_index;

        /**
         * @brief Shape index, e.g. the shape ID in shapegroup (if applicable).
         */
        uint32_t shape_index;

        /**
         * @brief Pointer to the associated shape.
         */
        const Shape* shape = nullptr;

        // Stores a pointer to the parent instance (if applicable)
        // Shape* instance = nullptr;

        /**
         * @brief Is the current interaction valid?
         * @return
         */
        bool isValid() const;

        /**
         * @brief Compute and return detailed information related to a surface interaction.
         * @param ray Ray to test against.
         * @param flags Information about what to query.
         * @return Surface interaction result.
         */
        SurfaceInteraction computeSurfaceInteraction(const Ray3d& ray, EHitComputeFlag flags);
    };
}
