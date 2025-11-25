#include <vislab/graphics/ray_differential.hpp>

namespace vislab
{
    RayDifferential3d::RayDifferential3d()
        : Ray3d()
        , origin_x(Eigen::Vector3d::Zero())
        , origin_y(Eigen::Vector3d::Zero())
        , direction_x(Eigen::Vector3d::Zero())
        , direction_y(Eigen::Vector3d::Zero())
        , hasDifferentials(false)
    {
    }

    RayDifferential3d::RayDifferential3d(const Ray3d& ray)
        : Ray3d(ray)
        , origin_x(Eigen::Vector3d::Zero())
        , origin_y(Eigen::Vector3d::Zero())
        , direction_x(Eigen::Vector3d::Zero())
        , direction_y(Eigen::Vector3d::Zero())
        , hasDifferentials(false)
    {
    }

    RayDifferential3d::RayDifferential3d(const Eigen::Vector3d& o, const Eigen::Vector3d& d)
        : origin_x(Eigen::Vector3d::Zero())
        , origin_y(Eigen::Vector3d::Zero())
        , direction_x(Eigen::Vector3d::Zero())
        , direction_y(Eigen::Vector3d::Zero())
        , hasDifferentials(false)
    {
        origin    = o;
        direction = d;
    }

    void RayDifferential3d::scaleDifferential(double amount)
    {
        origin_x    = (origin_x - origin) * amount + origin;
        origin_y    = (origin_y - origin) * amount + origin;
        direction_x = (direction_x - direction) * amount + direction;
        direction_y = (direction_y - direction) * amount + direction;
    }
}
