#include <vislab/graphics/rectangle.hpp>

#include <vislab/graphics/math.hpp>
#include <vislab/graphics/warp.hpp>

namespace vislab
{
    Rectangle::Rectangle()
        : mInverseSurfaceArea(0)
    {
        // update the internal precomputations whenever the matrix changed.
        this->transform.matrixChanged += [this](Transform* sender, const Eigen::Matrix4d* args)
        {
            this->update();
        };
        update();
    }

    void Rectangle::update()
    {
        Eigen::Vector3d dp_du  = transform.transformVector(Eigen::Vector3d(2, 0, 0));
        Eigen::Vector3d dp_dv  = transform.transformVector(Eigen::Vector3d(0, 2, 0));
        Eigen::Vector3d normal = transform.transformNormal(Eigen::Vector3d(0, 0, 1));
        normal.normalize();
        mFrame              = Frame(dp_du, dp_dv, normal);
        mInverseSurfaceArea = 1. / surfaceArea();
    }

    Eigen::AlignedBox3d Rectangle::objectBounds() const
    {
        return Eigen::AlignedBox3d(
            Eigen::Vector3d(-1, -1, 0),
            Eigen::Vector3d(1, 1, 0));
    }

    double Rectangle::surfaceArea() const
    {
        return mFrame.s.cross(mFrame.t).norm();
    }

    PositionSample Rectangle::samplePosition(const Eigen::Vector2d& sample) const
    {
        PositionSample ps;
        ps.position = transform.transformPoint(
            Eigen::Vector3d(sample.x() * 2.f - 1.f, sample.y() * 2.f - 1.f, 0.f));
        ps.normal = mFrame.n;
        ps.pdf    = mInverseSurfaceArea;
        ps.uv     = sample;
        ps.delta  = false;
        return ps;
    }

    double Rectangle::pdfPosition(const PositionSample& ps) const
    {
        return mInverseSurfaceArea;
    }

    PreliminaryIntersection Rectangle::preliminaryHit(const Ray3d& ray_) const
    {
        Ray ray               = transform.transformRayInverse(ray_);
        double t              = -ray.origin.z() / ray.direction.z();
        Eigen::Vector3d local = ray(t);

        // Is intersection within ray segment and rectangle?
        bool active = t >= ray.tMin && t <= ray.tMax && std::abs(local.x()) <= 1. && std::abs(local.y()) <= 1.;

        PreliminaryIntersection pi;
        pi.t          = active ? t : std::numeric_limits<double>::infinity();
        pi.prim_index = -1;
        pi.prim_uv    = local.xy();
        pi.shape      = this;
        return pi;
    }

    bool Rectangle::anyHit(const Ray3d& ray_) const
    {
        Ray3d ray               = transform.transformRayInverse(ray_);
        double t              = -ray.origin.z() / ray.direction.z();
        Eigen::Vector3d local = ray(t);

        // Is intersection within ray segment and rectangle?
        return t >= ray.tMin && t <= ray.tMax && std::abs(local.x()) <= 1. && std::abs(local.y()) <= 1.;
    }

    SurfaceInteraction Rectangle::computeSurfaceInteraction(const Ray3d& ray, const PreliminaryIntersection& pi, EHitComputeFlag flags) const
    {
        SurfaceInteraction si;
        Eigen::Vector2d prim_uv = pi.prim_uv;

        si.t = pi.t;
        // Re-project onto the rectangle to improve accuracy
        Eigen::Vector3d p = ray(pi.t);
        double dist       = (transform.getMatrix().translation() - p).dot(mFrame.n);
        si.position       = p + dist * mFrame.n;

        si.normal     = mFrame.n;
        si.sh_frame.n = mFrame.n;
        si.dp_du      = mFrame.s;
        si.dp_dv      = mFrame.t;
        si.uv         = Eigen::Vector2d(
            prim_uv.x() * 0.5f + 0.5f,
            prim_uv.y() * 0.5f + 0.5f);

        // si.dn_du = si.dn_dv = dr::zeros<Vector3f>();
        si.shape = this;

        // si.boundary_test = dr::min(0.5f - dr::abs(si.uv - 0.5f));

        return si;
    }
}
