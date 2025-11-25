#pragma once

#include <Eigen/Eigen>

namespace vislab
{
    class SurfaceInteraction;

    /**
     * @brief Position sample.
     */
    class PositionSample
    {
    public:
        /**
         * @brief Empty constructor.
         */
        PositionSample();

        /**
         * @brief Constructor.
         * @param si Surface interaction to extract position sample from.
         */
        PositionSample(const SurfaceInteraction& si);

        /**
         * @brief Constructor.
         * @param p Sampled position.
         * @param n Sampled surface normal (if applicable).
         * @param uv Optional UV coordiantes of the sampled position.
         * @param pdf Probability density at the sample.
         * @param delta Flag to determine if the sample was draw from a Dirac delta distribution.
         */
        PositionSample(const Eigen::Vector3d& p, const Eigen::Vector3d& n, const Eigen::Vector2d& uv, const double& pdf, const bool& delta);

        /**
         * @brief Sampled position.
         */
        Eigen::Vector3d position;

        /**
         * @brief Sampled surface normal (if applicable).
         */
        Eigen::Vector3d normal;

        /**
         * @brief Optional UV coordinates of the sampled position.
         */
        Eigen::Vector2d uv;

        /**
         * @brief Probability density at the sample.
         */
        double pdf;

        /**
         * @brief Flag to determine if the sample was drawn from a Dirac delta distribution.
         */
        bool delta;
    };
}
