#pragma once

#include "ray_fwd.hpp"

namespace vislab
{
    class Medium;

    /**
     * @brief Class that represents a ray.
     */
    template <typename TVector>
    class Ray
    {
    public:
        /**
         * @brief Vector type of the ray.
         */
        using Vector = TVector;

        /**
         * @brief Constructor.
         */
        Ray()
            : origin(Eigen::Vector3d::Zero())
            , direction(Eigen::Vector3d::Zero())
            , tMin(0)
            , tMax(std::numeric_limits<double>::max())
            , medium(nullptr)
        {
        }

        /**
         * @brief Constructor.
         * @param origin Origin of the ray.
         * @param direction Direction of the ray.
         * @param tMax Maximal length of the ray.
         * @param medium Medium that the ray starts in.
         */
        Ray(const Vector& origin, const Vector& direction, double tMin = 0, double tMax = std::numeric_limits<double>::max(), const Medium* medium = nullptr)
            : origin(origin)
            , direction(direction)
            , tMin(tMin)
            , tMax(tMax)
            , medium(medium)
        {
        }

        /**
         * @brief Samples the ray for a given distance using origin + t * direction.
         * @param t Distance on the ray.
         * @return Location on the ray.
         */
        Vector operator()(const double& t) const
        {
            return origin + t * direction;
        }

        /**
         * @brief Origin of the ray.
         */
        Vector origin;

        /**
         * @brief Direction of the ray.
         */
        Vector direction;

        /**
         * @brief Maximal length of the day.
         */
        mutable double tMin;

        /**
         * @brief Maximal length of the day.
         */
        mutable double tMax;

        /**
         * @brief Medium that the ray starts in.
         */
        const Medium* medium;
    };
}
