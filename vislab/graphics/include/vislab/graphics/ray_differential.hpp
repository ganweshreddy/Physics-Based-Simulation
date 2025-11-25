#pragma once

#include "ray.hpp"

namespace vislab
{
    /**
     * @brief Class the extends the ray with information about the neighboring pixels.
     */
    class RayDifferential3d : public Ray3d
    {
    public:
        /**
         * @brief Constructor.
         */
        RayDifferential3d();

        /**
         * @brief Construct from a Ray instance
         * @param ray Ray to initalize from.
         */
        RayDifferential3d(const Ray3d& ray);

        /**
         * @brief Construct a new ray (o, d)
         * @param o Origin.
         * @param d Direction
         */
        RayDifferential3d(const Eigen::Vector3d& o, const Eigen::Vector3d& d);

        /**
         * @brief Scales the differentials by a certain amount.
         * @param amount Amount to scale differentials by.
         */
        void scaleDifferential(double amount);

        /**
         * @brief Origin of the right neighbor.
         */
        Eigen::Vector3d origin_x;

        /**
         * @brief Origin of the top neighbor.
         */
        Eigen::Vector3d origin_y;

        /**
         * @brief Direction of the ray of the right neighbor.
         */
        Eigen::Vector3d direction_x;

        /**
         * @brief Direction of the ray of the top neighbor.
         */
        Eigen::Vector3d direction_y;

        /**
         * @brief Flag that determines whether differentials are defined yet.
         */
        bool hasDifferentials;
    };
}
