#include <vislab/graphics/sphere.hpp>

#include <vislab/graphics/math.hpp>
#include <vislab/graphics/warp.hpp>

namespace vislab
{
    Sphere::Sphere()
        : mInverseSurfaceArea(0)
    {
        // update the internal precomputations whenever the matrix changed.
        this->transform.matrixChanged += [this](Transform* sender, const Eigen::Matrix4d* args)
        {
            this->update();
        };
        update();
    }

    void Sphere::update()
    {
        mInverseSurfaceArea = 1. / surfaceArea();
    }

    Eigen::AlignedBox3d Sphere::objectBounds() const
    {
        return Eigen::AlignedBox3d(
            Eigen::Vector3d(-1, -1, -1),
            Eigen::Vector3d(1, 1, 1));
    }

    double Sphere::surfaceArea() const
    {
        return 4.f * EIGEN_PI * radius() * radius();
    }

    double Sphere::radius() const
    {
        return this->transform.getMatrix()(0, 0);
    }

    Eigen::Vector3d Sphere::center() const
    {
        return this->transform.getMatrix().translation();
    }

    PositionSample Sphere::samplePosition(const Eigen::Vector2d& sample) const
    {
        Eigen::Vector3d local = Warp::square_to_uniform_sphere(sample);

        PositionSample ps;
        ps.position = local * radius() + center();
        ps.normal   = local;

        if (reverseOrientation)
            ps.normal = -ps.normal;

        ps.delta = radius() == 0.f;
        ps.pdf   = mInverseSurfaceArea;
        return ps;
    }

    double Sphere::pdfPosition(const PositionSample& ps) const
    {
        return mInverseSurfaceArea;
    }

    DirectionSample Sphere::sampleDirection(const Interaction& it, const Eigen::Vector2d& sample) const
    {
        DirectionSample result;

        Eigen::Vector3d dc_v = center() - it.position;
        double dc_2          = dc_v.squaredNorm();

        double radius_adj = radius() * (reverseOrientation ? (1. + std::numeric_limits<double>::epsilon()) : (1.f - std::numeric_limits<double>::epsilon()));
        bool outside_mask = dc_2 > radius_adj * radius_adj;
        if (outside_mask)
        {
            double inv_dc            = 1 / std::sqrt(dc_2),
                   sin_theta_max     = radius() * inv_dc,
                   sin_theta_max_2   = sin_theta_max * sin_theta_max,
                   inv_sin_theta_max = 1.0 / sin_theta_max,
                   cos_theta_max     = std::sqrt(std::max(0., 1. - sin_theta_max_2));

            /* Fall back to a Taylor series expansion for small angles, where
               the standard approach suffers from severe cancellation errors */
            double sin_theta_2 = sin_theta_max_2 > 0.00068523f ? /* sin^2(1.5 deg) */
                                     1. - ((cos_theta_max - 1.) * sample.x() + 1.) * ((cos_theta_max - 1.) * sample.x() + 1.)
                                                               : sin_theta_max_2 * sample.x(),
                   cos_theta   = std::sqrt(std::max(0., 1. - sin_theta_2));

            // Based on https://www.akalin.com/sampling-visible-sphere
            double cos_alpha = sin_theta_2 * inv_sin_theta_max +
                               cos_theta * std::sqrt(std::max(0., 1. - sin_theta_2 * (inv_sin_theta_max * inv_sin_theta_max))),
                   sin_alpha = std::sqrt(std::max(0., 1 - cos_alpha * cos_alpha));

            double sin_phi = std::sin(sample.y() * (2.f * EIGEN_PI)),
                   cos_phi = std::cos(sample.y() * (2.f * EIGEN_PI));

            Eigen::Vector3d d = Frame(dc_v * -inv_dc).toWorld(Eigen::Vector3d(cos_phi * sin_alpha, sin_phi * sin_alpha, cos_alpha));

            DirectionSample ds;
            ds.position  = d * radius() + center();
            ds.normal    = d;
            ds.direction = ds.position - it.position;

            double dist2 = ds.direction.squaredNorm();
            ds.dist      = sqrt(dist2);
            ds.direction = ds.direction / ds.dist;
            ds.pdf       = Warp::square_to_uniform_cone_pdf(Eigen::Vector3d(0, 0, 0), cos_theta_max);
            if (ds.dist == 0.)
                ds.pdf = 0;

            result = ds;
        }

        bool inside_mask = !outside_mask;
        if (inside_mask)
        {
            Eigen::Vector3d d = Warp::square_to_uniform_sphere(sample);
            DirectionSample ds;
            ds.position  = d * radius() + center();
            ds.normal    = d;
            ds.direction = ds.position - it.position;

            double dist2 = ds.direction.squaredNorm();
            ds.dist      = sqrt(dist2);
            ds.direction = ds.direction / ds.dist;
            ds.pdf       = mInverseSurfaceArea * dist2 / std::abs(ds.direction.dot(ds.normal));

            result = ds;
        }

        result.delta = radius() == 0.f;

        if (reverseOrientation)
            result.normal = -result.normal;

        return result;
    }

    double Sphere::pdfDirection(const Interaction& it, const DirectionSample& ds) const
    {
        // Sine of the angle of the cone containing the sphere as seen from 'it.p'.
        double sin_alpha = radius() * 1. / (center() - it.position).norm(),
               cos_alpha = std::sqrt(std::max(0., 1.f - sin_alpha * sin_alpha));

        return sin_alpha < 1 - std::numeric_limits<double>::epsilon() ?
                                                                      // Reference point lies outside the sphere
                   Warp::square_to_uniform_cone_pdf(Eigen::Vector3d(0, 0, 0), cos_alpha)
                                                                      : mInverseSurfaceArea * (ds.dist * ds.dist) / std::abs(ds.direction.dot(ds.normal));
    }

    PreliminaryIntersection Sphere::preliminaryHit(const Ray3d& ray) const
    {
        double mint = ray.tMin;
        double maxt = ray.tMax;

        Eigen::Vector3d o = ray.origin - center();
        Eigen::Vector3d d(ray.direction);

        double A = d.squaredNorm();
        double B = 2. * o.dot(d);
        double C = o.squaredNorm() - radius() * radius();

        auto [solution_found, near_t, far_t] = solveQuadratic(A, B, C);

        // Sphere doesn't intersect with the segment on the ray
        bool out_bounds = !(near_t <= maxt && far_t >= mint); // NaN-aware conditionals

        // Sphere fully contains the segment of the ray
        bool in_bounds = near_t < mint && far_t > maxt;

        bool active = solution_found && !out_bounds && !in_bounds;

        PreliminaryIntersection pi;
        pi.t     = active ? (near_t < mint ? far_t : near_t) : std::numeric_limits<double>::infinity();
        pi.shape = this;

        return pi;
    }

    bool Sphere::anyHit(const Ray3d& ray) const
    {
        double mint = ray.tMin;
        double maxt = ray.tMax;

        Eigen::Vector3d o = ray.origin - center();
        Eigen::Vector3d d(ray.direction);

        double A = d.squaredNorm();
        double B = 2. * o.dot(d);
        double C = o.squaredNorm() - radius() * radius();

        auto [solution_found, near_t, far_t] = solveQuadratic(A, B, C);

        // Sphere doesn't intersect with the segment on the ray
        bool out_bounds = !(near_t <= maxt && far_t >= mint); // NaN-aware conditionals

        // Sphere fully contains the segment of the ray
        bool in_bounds = near_t < mint && far_t > maxt;

        return solution_found && !out_bounds && !in_bounds;
    }

    SurfaceInteraction Sphere::computeSurfaceInteraction(const Ray3d& ray, const PreliminaryIntersection& pi, EHitComputeFlag flags) const
    {
        bool active = pi.isValid();

        SurfaceInteraction si;
        si.t = active ? pi.t : std::numeric_limits<double>::infinity();

        si.sh_frame.n = (ray(pi.t) - center()).normalized();

        // Re-project onto the sphere to improve accuracy
        si.position = si.sh_frame.n * radius() + center();

        if (hasFlag(flags, EHitComputeFlag::UV))
        {
            Eigen::Vector3d local = (transform.getMatrixInverse() * Eigen::Vector4d(si.position.x(), si.position.y(), si.position.z(), 1)).xyz();

            double rd_2  = local.x() * local.x() + local.y() * local.y(),
                   theta = unitAngleZ(local),
                   phi   = atan2(local.y(), local.x());

            if (phi < 0.)
                phi += 2.f * EIGEN_PI;

            si.uv = Eigen::Vector2d(phi * 0.5 / EIGEN_PI, theta / EIGEN_PI);

            if (hasFlag(flags, EHitComputeFlag::dPdUV))
            {
                si.dp_du = Eigen::Vector3d(-local.y(), local.x(), 0.f);

                double rd      = std::sqrt(rd_2),
                       inv_rd  = 1. / rd,
                       cos_phi = local.x() * inv_rd,
                       sin_phi = local.y() * inv_rd;

                si.dp_dv = Eigen::Vector3d(local.z() * cos_phi,
                                           local.z() * sin_phi,
                                           -rd);

                bool singularity_mask = active && rd == 0.;
                if (singularity_mask)
                    si.dp_dv = Eigen::Vector3d(1.f, 0.f, 0.f);

                si.dp_du = transform.getMatrix().block(0, 0, 3, 3) * si.dp_du * (2.f * EIGEN_PI);
                si.dp_dv = transform.getMatrix().block(0, 0, 3, 3) * si.dp_dv * EIGEN_PI;
            }
        }

        if (reverseOrientation)
            si.sh_frame.n = -si.sh_frame.n;

        si.normal = si.sh_frame.n;

        return si;
    }
}
