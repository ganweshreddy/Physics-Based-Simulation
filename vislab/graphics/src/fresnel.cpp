#include <vislab/graphics/fresnel.hpp>

#include <vislab/graphics/math.hpp>

namespace vislab
{
    std::tuple<double, double, double, double> fresnel(double cos_theta_i, double eta)
    {
        auto outside_mask = cos_theta_i >= 0.f;

        double rcp_eta = 1. / eta,
               eta_it  = outside_mask ? eta : rcp_eta,
               eta_ti  = outside_mask ? rcp_eta : eta;

        // Using Snell's law, calculate the squared sine of the angle between the surface normal and the transmitted ray
        double cos_theta_t_sqr =
            (-(-cos_theta_i * cos_theta_i + 1.f) * eta_ti * eta_ti + 1.f);

        // Find the absolute cosines of the incident/transmitted rays
        double cos_theta_i_abs = abs(cos_theta_i);
        double cos_theta_t_abs = safe_sqrt(cos_theta_t_sqr);

        auto index_matched = eta == 1.,
             special_case  = index_matched || cos_theta_i_abs == 0.;

        double r_sc = index_matched ? 0. : 1.;

        // Amplitudes of reflected waves
        double a_s = (-eta_it * cos_theta_t_abs + cos_theta_i_abs) /
                     (eta_it * cos_theta_t_abs + cos_theta_i_abs);

        double a_p = (-eta_it * cos_theta_i_abs + cos_theta_t_abs) /
                     (eta_it * cos_theta_i_abs + cos_theta_t_abs);

        double r = .5 * ((a_s * a_s) + (a_p * a_p));

        if (special_case)
            r = r_sc;

        // Adjust the sign of the transmitted direction
        double cos_theta_t = cos_theta_t_abs * (cos_theta_i < 0 ? 1 : -1);

        return { r, cos_theta_t, eta_it, eta_ti };
    }

    Eigen::Vector3d reflect(const Eigen::Vector3d& wi)
    {
        return Eigen::Vector3d(-wi.x(), -wi.y(), wi.z());
    }

    Eigen::Vector3d reflect(const Eigen::Vector3d& wi, const Eigen::Vector3d& n)
    {
        return n * 2.f * wi.dot(n) - wi;
    }

    Eigen::Vector3d refract(const Eigen::Vector3d& wi, double cos_theta_t, double eta_ti)
    {
        return Eigen::Vector3d(-eta_ti * wi.x(), -eta_ti * wi.y(), cos_theta_t);
    }

    Eigen::Vector3d refract(const Eigen::Vector3d& wi, const Eigen::Vector3d& n, double cos_theta_t, double eta_ti)
    {
        return (n * (wi.dot(n) * eta_ti + cos_theta_t) - wi * eta_ti);
    }
}
