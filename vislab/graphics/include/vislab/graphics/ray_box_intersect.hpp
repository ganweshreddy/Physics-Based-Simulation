#pragma once

#include "ray_fwd.hpp"

namespace vislab
{
    /**
     * @brief Ray box intersection test.
     */
    class RayBoxIntersection
    {
    public:
        /**
         * @brief Tests if a ray has any intersection with an axis aligned bounding box.
         * @param ray Ray to test.
         * @param bounds Bounding box to test against.
         * @return True if the ray intersects the box.
         */
        static bool anyHit(const Ray2d& ray, const Eigen::AlignedBox2d& bounds);

        /**
         * @brief Tests if a ray has any intersection with an axis aligned bounding box.
         * @param ray Ray to test.
         * @param bounds Bounding box to test against.
         * @return True if the ray intersects the box.
         */
        static bool anyHit(const Ray3d& ray, const Eigen::AlignedBox3d& bounds);

        /**
         * @brief Computes the intersection of a ray with an axis aligned bounding box.
         * @param ray Ray to test.
         * @param bounds Bounding box to test against.
         * @param t_out Entry and exit distance.
         * @return True if the ray intersects the box.
         */
        static bool closestHit(const Ray2d& ray, const Eigen::AlignedBox2d& bounds, Eigen::Vector2d& t_out);

        /**
         * @brief Computes the intersection of a ray with an axis aligned bounding box.
         * @param ray Ray to test.
         * @param bounds Bounding box to test against.
         * @param t_out Entry and exit distance.
         * @return True if the ray intersects the box.
         */
        static bool closestHit(const Ray3d& ray, const Eigen::AlignedBox3d& bounds, Eigen::Vector2d& t_out);
    };
}
