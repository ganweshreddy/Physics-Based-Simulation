#include <vislab/graphics/ray_line_intersect.hpp>

#include <vislab/graphics/ray.hpp>

namespace vislab
{
    bool RayLineIntersection::anyHit(const Ray2d& ray, const Eigen::Vector2d& point1, const Eigen::Vector2d& point2)
    {
        Eigen::Vector2d v1 = ray.origin - point1;
        Eigen::Vector2d v2 = point2 - point1;
        Eigen::Vector2d v3 = Eigen::Vector2d(-ray.direction.y(), ray.direction.x());

        double dot = v2.dot(v3);
        if (std::abs(dot) < std::numeric_limits<double>::epsilon())
            return false;

        double t2 = v1.dot(v3) / dot;
        if (t2 >= 0.0 && t2 <= 1.0)
        {
            double t1 = (v2.x() * v1.y() - v2.y() * v1.x()) / dot;
            return ray.tMin < t1 && t1 < ray.tMax;
        }
        return false;
    }

    bool RayLineIntersection::closestHit(const Ray2d& ray, const Eigen::Vector2d& point1, const Eigen::Vector2d& point2, double& t_out)
    {
        Eigen::Vector2d v1 = ray.origin - point1;
        Eigen::Vector2d v2 = point2 - point1;
        Eigen::Vector2d v3 = Eigen::Vector2d(-ray.direction.y(), ray.direction.x());

        double dot = v2.dot(v3);
        if (std::abs(dot) < std::numeric_limits<double>::epsilon())
            return false;

        double t2 = v1.dot(v3) / dot;
        if (t2 >= 0.0 && t2 <= 1.0)
        {
            double t1 = (v2.x() * v1.y() - v2.y() * v1.x()) / dot;
            if (ray.tMin < t1 && t1 < ray.tMax)
            {
                t_out = t1;
                return true;
            }
        }
        return false;
    }
}
