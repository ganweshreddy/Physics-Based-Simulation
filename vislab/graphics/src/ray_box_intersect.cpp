#include <vislab/graphics/ray_box_intersect.hpp>

#include <vislab/graphics/ray.hpp>

namespace vislab
{
    static bool ray_check(const Eigen::Vector2d& t, double tmin, double tmax)
    {
        return !((t.x() <= tmin && t.y() <= tmin) || (tmax <= t.x() && tmax <= t.y()));
    }

    bool RayBoxIntersection::anyHit(const Ray2d& ray, const Eigen::AlignedBox2d& bounds)
    {
        Eigen::Vector2d t1 = (bounds.min() - ray.origin).cwiseQuotient(ray.direction);
        Eigen::Vector2d t2 = (bounds.max() - ray.origin).cwiseQuotient(ray.direction);
        double tmin        = std::max(std::min(t1.x(), t2.x()), std::min(t1.y(), t2.y()));
        double tmax        = std::min(std::max(t1.x(), t2.x()), std::max(t1.y(), t2.y()));
        return tmin <= tmax && ray_check(Eigen::Vector2d(tmin, tmax), ray.tMin, ray.tMax);
    }
    
    bool RayBoxIntersection::anyHit(const Ray3d& ray, const Eigen::AlignedBox3d& bounds)
    {
        Eigen::Vector2d t;
        return closestHit(ray, bounds, t);
    }

    bool RayBoxIntersection::closestHit(const Ray2d& ray, const Eigen::AlignedBox2d& bounds, Eigen::Vector2d& t_out)
    {
        Eigen::Vector2d t1 = (bounds.min() - ray.origin).cwiseQuotient(ray.direction);
        Eigen::Vector2d t2 = (bounds.max() - ray.origin).cwiseQuotient(ray.direction);
        double tmin        = std::max(std::min(t1.x(), t2.x()), std::min(t1.y(), t2.y()));
        double tmax        = std::min(std::max(t1.x(), t2.x()), std::max(t1.y(), t2.y()));
        if (tmin <= tmax)
        {
            t_out = Eigen::Vector2d(tmin, tmax);
            return ray_check(t_out, ray.tMin, ray.tMax);
        }
        return false;
    }

    bool RayBoxIntersection::closestHit(const Ray3d& ray, const Eigen::AlignedBox3d& bounds, Eigen::Vector2d& t_out)
    {
        Eigen::Vector3d tmin = (bounds.min() - ray.origin).cwiseQuotient(ray.direction);
        Eigen::Vector3d tmax = (bounds.max() - ray.origin).cwiseQuotient(ray.direction);

        if (tmin.x() > tmax.x())
            std::swap(tmin.x(), tmax.x());
        if (tmin.y() > tmax.y())
            std::swap(tmin.y(), tmax.y());
        t_out = Eigen::Vector2d(tmin.x(), tmax.x());

        if ((t_out.x() > tmax.y()) || (tmin.y() > t_out.y()))
            return false;

        if (tmin.y() > t_out.x())
            t_out.x() = tmin.y();
        if (tmax.y() < t_out.y())
            t_out.y() = tmax.y();

        if (tmin.z() > tmax.z())
            std::swap(tmin.z(), tmax.z());
        if ((tmin.x() > tmax.z()) || (tmin.z() > tmax.x()))
            return false;

        if (tmin.z() > t_out.x())
            t_out.x() = tmin.z();
        if (tmax.z() < t_out.y())
            t_out.y() = tmax.z();

        return ray_check(t_out, ray.tMin, ray.tMax);
    }
}
