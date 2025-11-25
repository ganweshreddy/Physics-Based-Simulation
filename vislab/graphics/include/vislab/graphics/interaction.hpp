#pragma once

#include <Eigen/Eigen>

namespace vislab
{
    /**
     * @brief Represents the interaction along a ray.
     */
    class Interaction
    {
    public:
        /**
         * @brief Constructor.
         */
        Interaction();

        /**
         * @brief Constructor.
         * @param tt Traveled distance along the ray.
         * @param pos Position of the interaction in world coordinates.
         */
        Interaction(double t, const Eigen::Vector3d& pos);

        /**
         * @brief Is the current interaction valid?
         * @return True if valid.
         */
        bool isValid() const;

        /**
         * @brief Distance traveled along the ray.
         */
        double t;

        /**
         * @brief Position of the interaction in world coordinates.
         */
        Eigen::Vector3d position;
    };
}
