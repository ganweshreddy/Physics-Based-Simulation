#include <vislab/graphics/shape.hpp>

#include <vislab/graphics/bsdf.hpp>

namespace vislab
{
    Shape::Shape()
        : reverseOrientation(false)
    {
    }

    Eigen::AlignedBox3d Shape::worldBounds() const
    {
        Eigen::AlignedBox3d oob = objectBounds();
        Eigen::AlignedBox3d aabb;
        aabb.setEmpty();
        aabb.extend(transform.transformPoint(oob.corner(Eigen::AlignedBox3d::CornerType::BottomLeftFloor)));
        aabb.extend(transform.transformPoint(oob.corner(Eigen::AlignedBox3d::CornerType::BottomRightFloor)));
        aabb.extend(transform.transformPoint(oob.corner(Eigen::AlignedBox3d::CornerType::TopLeftFloor)));
        aabb.extend(transform.transformPoint(oob.corner(Eigen::AlignedBox3d::CornerType::TopRightFloor)));
        aabb.extend(transform.transformPoint(oob.corner(Eigen::AlignedBox3d::CornerType::BottomLeftCeil)));
        aabb.extend(transform.transformPoint(oob.corner(Eigen::AlignedBox3d::CornerType::BottomRightCeil)));
        aabb.extend(transform.transformPoint(oob.corner(Eigen::AlignedBox3d::CornerType::TopLeftCeil)));
        aabb.extend(transform.transformPoint(oob.corner(Eigen::AlignedBox3d::CornerType::TopRightCeil)));
        return aabb;
    }

    DirectionSample Shape::sampleDirection(const Interaction& it, const Eigen::Vector2d& sample) const
    {
        DirectionSample ds(samplePosition(sample));
        ds.direction = ds.position - it.position;

        double dist_squared = ds.direction.squaredNorm();
        ds.dist             = sqrt(dist_squared);
        ds.direction /= ds.dist;

        double dp = std::abs(ds.direction.dot(ds.normal));
        ds.pdf *= dp == 0. ? 0. : dist_squared / dp;
        // ds.object = (const Object*)this;

        return ds;
    }

    double Shape::pdfDirection(const Interaction& it, const DirectionSample& ds) const
    {
        double pdf = pdfPosition(ds),
               dp  = std::abs(ds.direction.dot(ds.normal));

        pdf *= dp == 0. ? 0. : (ds.dist * ds.dist) / dp;
        return pdf;
    }

    bool Shape::anyHit(const Ray3d& ray) const
    {
        return preliminaryHit(ray).isValid();
    }

    int32_t Shape::countHits(const Ray3d& _ray) const
    {
        Ray3d ray             = _ray;
        int32_t intersections = 0;
        while (true)
        {
            auto pi = preliminaryHit(ray);
            if (!pi.isValid())
                break;
            ray.tMin = pi.t + rayEpsilon();
            intersections++;
        }
        return intersections;
    }

    SurfaceInteraction Shape::closestHit(const Ray3d& ray, EHitComputeFlag flags) const
    {
        auto pi = preliminaryHit(ray);
        return pi.computeSurfaceInteraction(ray, flags);
    }

    void Shape::buildAccelerationTree()
    {
    }
}
