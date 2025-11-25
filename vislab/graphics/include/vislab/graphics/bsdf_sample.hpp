#pragma once

#include "bsdf_flag.hpp"

#include <Eigen/Eigen>

namespace vislab
{
    /**
     * @brief Data structure holding the result of BSDF sampling operations.
     */
    struct BSDFSample
    {
        /**
         * @brief Normalized outgoing direction in local coordinates.
         */
        Eigen::Vector3d wo;

        /**
         * @brief Probability density at the sample.
         */
        double pdf;

        /**
         * @brief Relative index of refraction in the sampled direction.
         */
        double eta;

        /**
         * @brief Stores the component type that was sampled by \ref BSDF::sample()
         */
        uint32_t sampledType;

        /**
         * @brief Tests if a flag is set.
         * @param flag Flag to test.
         * @return True, if the flag is set.
         */
        bool hasFlag(const EBSDFFlag& flag) const;

        /**
         * @brief Constructor.
         * @param wo Outgoing direction in local coordinates.
         */
        BSDFSample(const Eigen::Vector3d& wo);
    };
}
