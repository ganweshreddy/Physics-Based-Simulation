#pragma once

#include "position_sample.hpp"

namespace vislab
{
    class Interaction;
    class Light;
    class Scene;

    /**
     * @brief Direction sample.
     */
    class DirectionSample : public PositionSample
    {
    public:
        /**
         * @brief Empty constructor.
         */
        DirectionSample();

        /**
         * @brief Constructor.
         * @param si Surface interaction at which direction sample starts.
         * @param ref Interaction at which direction sample ends.
         */
        DirectionSample(const Scene* scene, const SurfaceInteraction& si, const Interaction& ref);

        /**
         * @brief Constructor.
         * @param p Sampled position.
         * @param n Sampled surface normal (if applicable).
         * @param uv Optional UV coordiantes of the sampled position.
         * @param pdf Probability density at the sample.
         * @param delta Flag to determine if the sample was draw from a Dirac delta distribution.
         * @param d Unit direction from the reference point to the target shape.
         * @param dist Distance from the reference point to the target shape.
         */
        DirectionSample(const Eigen::Vector3d& p, const Eigen::Vector3d& n, const Eigen::Vector2d& uv, const double& pdf, const bool& delta, const Eigen::Vector3d& d, const double& dist, const Light* light);

        /**
         * @brief Constructor.
         * @param base Base class to build object from.
         */
        DirectionSample(const PositionSample& base);

        /**
         * @brief Unit direction from the reference point to the target shape.
         */
        Eigen::Vector3d direction;

        /**
         * @brief Distance from the reference point to the target shape.
         */
        double dist;

        /**
         * @brief Optional: pointer to a light source if this direction was sampled to a light source.
         */
        const Light* light;
    };
}
