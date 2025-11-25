#pragma once

#include <Eigen/Eigen>

#include <tuple>

namespace vislab
{
    /**
     * @brief Unpolarized Fresnel reflection coefficient
     * @param cos_theta_i Cosine of the angle between the surface normal and the incident ray.
     * @param eta Relative refractive index of the interface.
     * @return A tuple consisting of Fresnel reflection coefficient, cosine of the angle between the surface normal and the transmitted ray, relative index of refraction in the direction of travel, reciprocal of the relative index of refraction in the direction of travel.
     */
    std::tuple<double, double, double, double> fresnel(double cos_theta_i, double eta);

    /**
     * @brief Reflection in local coordinates
     * @param wi Indicident direction.
     * @return Reflected vector.
     */
    Eigen::Vector3d reflect(const Eigen::Vector3d& wi);

    /**
     * @brief Reflect wi with respect to a given surface normal n.
     * @param wi Incident direction.
     * @param n Surface normal.
     * @return Reflected vector.
     */
    Eigen::Vector3d reflect(const Eigen::Vector3d& wi, const Eigen::Vector3d& n);

    /**
     * @brief Refraction in local coordinates.
     * @param wi Incident direction.
     * @param cos_theta_t Cosine of the angle between the normal the transmitted ray.
     * @param eta_ti Relative index of refraction (transmitted / incident).
     * @return Refracted vector.
     */
    Eigen::Vector3d refract(const Eigen::Vector3d& wi, double cos_theta_t, double eta_ti);

    /**
     * @brief Refract wi with respect to a given surface normal.
     * @param wi Incident direction.
     * @param n Surface normal.
     * @param cos_theta_t Cosine of the angle between the normal the transmitted ray.
     * @param eta_ti Relative index of refraction (transmitted / incident).
     * @return Refracted vector.
     */
    Eigen::Vector3d refract(const Eigen::Vector3d& wi, const Eigen::Vector3d& n, double cos_theta_t, double eta_ti);
}
