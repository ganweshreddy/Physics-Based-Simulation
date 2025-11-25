#pragma once

#include "math.hpp"

#include <Eigen/Eigen>

namespace vislab
{
    /**
     * @brief Herlper functions that warp between different spaces.
     */
    class Warp
    {
    public:
        /**
         * @brief Transformation from a unit square [0,1]^2 to a uniform sphere.
         * @param sample Unit square sample.
         * @return Uniform sphere sample.
         */
        static Eigen::Vector3d square_to_uniform_sphere(const Eigen::Vector2d& sample)
        {
            double z    = -2. * sample.y() + 1.,
                   r    = circ(z);
            auto [s, c] = sincos(2. * EIGEN_PI * sample.x());
            return Eigen::Vector3d(r * c, r * s, z);
        }

        /**
         * @brief Uniformly sample a vector that lies within a given cone of angles around the Z axis.
         * @param sample Unit square sample.
         * @param cos_cutoff Cutoff for cosine.
         * @return Uniform vector in cone.
         */
        static Eigen::Vector3d square_to_uniform_cone(const Eigen::Vector2d& sample, const double& cos_cutoff)
        {
#if 0
			// Approach 1: warping method based on standard disk mapping
			double cos_theta = (1. - sample.y()) + sample.y() * cos_cutoff,
				sin_theta = circ(cos_theta);

			auto [s, c] = sincos(2 * M_PI * sample.x());
			return Eigen::Vector3d(c * sin_theta, s * sin_theta, cos_theta);
#else
            // Approach 2: low-distortion warping technique based on concentric disk mapping
            double one_minus_cos_cutoff(1. - cos_cutoff);
            Eigen::Vector2d p = square_to_uniform_disk_concentric(sample);
            double pn         = p.squaredNorm();
            double z          = cos_cutoff + one_minus_cos_cutoff * (1.f - pn);
            p *= safe_sqrt(one_minus_cos_cutoff * (2.f - one_minus_cos_cutoff * pn));
            return Eigen::Vector3d(p.x(), p.y(), z);
#endif
        }

        /**
         * @brief Density of square_to_uniform_cone per unit area.
         * @param v Position to compute density for.
         * @param cos_cutoff Cosine cutoff.
         * @return Density of sample.
         */
        static double square_to_uniform_cone_pdf(const Eigen::Vector3d& v, const double& cos_cutoff)
        {
            return 0.5 / EIGEN_PI / (1. - cos_cutoff);
        }

        /**
         * @brief Low-distortion concentric square to disk mapping by Peter Shirley based on the concentric map code with less branching by Dave Cline, see http://psgraphics.blogspot.ch/2011/01/improved-code-for-concentric-map.html
         * @param sample Uniform square sample.
         * @return Disk sample.
         */
        static Eigen::Vector2d square_to_uniform_disk_concentric(const Eigen::Vector2d& sample)
        {
            double x = 2. * sample.x() - 1.,
                   y = 2. * sample.y() - 1.;

            bool is_zero         = x == 0 && y == 0,
                 quadrant_1_or_3 = abs(x) < abs(y);

            double r  = quadrant_1_or_3 ? y : x,
                   rp = quadrant_1_or_3 ? x : y;

            double phi = .25 * EIGEN_PI * rp / r;
            if (quadrant_1_or_3)
                phi = .5 * EIGEN_PI - phi;
            if (is_zero)
                phi = 0;

            auto [s, c] = sincos(phi);
            return Eigen::Vector2d(r * c, r * s);
        }

        /**
         * @brief Sample a cosine-weighted vector on the unit hemisphere with respect to solid angles.
         * @param sample Sample point in [0,1]^2.
         * @return Sample point on hemisphere.
         */
        static Eigen::Vector3d square_to_cosine_hemisphere(const Eigen::Vector2d& sample)
        {
            // Low-distortion warping technique based on concentric disk mapping
            Eigen::Vector2d p = square_to_uniform_disk_concentric(sample);

            // Guard against numerical imprecisions
            double z = safe_sqrt(1. - p.squaredNorm());

            return Eigen::Vector3d(p.x(), p.y(), z);
        }

        /**
         * @brief Density of \ref square_to_cosine_hemisphere() with respect to solid angles
         * @tparam TestDomain Safety checks enabled.
         * @param v Sample point.
         * @return Probability density.
         */
        template <bool TestDomain = false>
        static double square_to_cosine_hemisphere_pdf(const Eigen::Vector3d& v)
        {
            const double math_ray_epsilon = std::numeric_limits<double>::epsilon();
            if constexpr (TestDomain)
                return (std::abs(v.squaredNorm() - 1.) > math_ray_epsilon ||
                        v.z() < 0.)
                           ? 0.
                           : (1.0 / EIGEN_PI * v.z());
            else
                return 1.0 / EIGEN_PI * v.z();
        }

    private:
        /**
         * @brief Circle arc.
         * @param x Value to sample y coordinate for.
         * @return y coordinate of circle arc.
         */
        static double circ(double x) { return safe_sqrt(1. - x * x); }

        /**
         * @brief Joint computation of sine and cosine.
         * @param x Value to compute sine and cosine for.
         * @return Sine and cosine.
         */
        static std::pair<double, double> sincos(double x) { return std::make_pair<double, double>(std::sin(x), std::cos(x)); }
    };
}