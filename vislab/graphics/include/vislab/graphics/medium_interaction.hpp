#pragma once

#include "interaction.hpp"

#include "frame.hpp"

namespace vislab
{
    class Medium;

    /**
     * @brief Represents an interaction in a medium.
     */
    class MediumInteraction : public Interaction
    {
    public:
        /**
         * @brief Medium that was interacted with.
         */
        const Medium* medium;

        /**
         * @brief Shading frame.
         */
        Frame sh_frame;

        /**
         * @brief Incident ray direction in the local shading frame.
         */
        Eigen::Vector3d wi;

        /**
         * @brief Convert a local shading-space vector into world space.
         * @param v Local shading-space vector to transform.
         * @return Transformed vector in world space.
         */
        Eigen::Vector3d toWorld(const Eigen::Vector3d& v) const;

        /**
         * @brief Convert a world-space vector into local shading coordinates.
         * @param v World-space vector to transform.
         * @return Transformed vector in local shading coordinates.
         */
        Eigen::Vector3d toLocal(const Eigen::Vector3d& v) const;
    };
}
