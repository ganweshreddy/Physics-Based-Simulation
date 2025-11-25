#pragma once

#include "ray_fwd.hpp"

namespace vislab
{
    /**
     * @brief Ray line intersection test.
     */
    class RayLineIntersection
    {
    public:
        /**
         * @brief Tests if a ray intersects a line segment defined by its end points point1 and point2.
         * @param ray Ray to test.
         * @param point1 Start point of the line segment.
         * @param point2 End point of the line segment.
         * @return True if an intersection occured.
         */
        static bool anyHit(const Ray2d& ray, const Eigen::Vector2d& point1, const Eigen::Vector2d& point2);

        /**
         * @brief Tests if a ray intersects a line segment defined by its end points point1 and point2.
         * @param ray Ray to test.
         * @param point1 Start point of the line segment.
         * @param point2 End point of the line segment.
         * @param t_out Ray parameter at which intersection occured.
         * @return True if an intersection occured.
         */
        static bool closestHit(const Ray2d& ray, const Eigen::Vector2d& point1, const Eigen::Vector2d& point2, double& t_out);
    };
}
