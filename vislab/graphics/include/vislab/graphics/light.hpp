#pragma once

#include "endpoint.hpp"
#include "light_flag.hpp"
#include "transform.hpp"

#include <Eigen/Eigen>

namespace vislab
{
    /**
     * @brief Interface for light sources.
     */
    class Light : public Endpoint
    {
        VISLAB_INTERFACE(Light, Data)

    public:
        /**
         * @brief Is the position a dirac delta?
         * @return True, if position is dirac delta.
         */
        bool isDeltaPosition() const;

        /**
         * @brief Is the direction a dirac delta?
         * @return True, if direction is dirac delta.
         */
        bool isDeltaDirection() const;

        /**
         * @brief Is an area light source?
         * @return True, if area light source.
         */
        bool isArea() const;

        /**
         * @brief Is an infinite light source?
         * @return True, if infinite light source.
         */
        bool isInfinite() const;

        /**
         * @brief Is this an environment map light emitter?
         * @return True, if an environment map light emitter.
         */
        bool isEnvironment() const;

        /**
         * @brief Is isDeltaPosition or isDeltaDirection?
         * @return True, if is delta.
         */
        bool isDeltaLight() const;

        /**
         * @brief Transformation of a shape from object space to world space.
         */
        Transform transform;

        /**
         * @brief Tests if a flag is set.
         * @param flag Flag to test.
         * @return True, if the flag is set.
         */
        bool hasFlag(const ELightFlag& flag) const;

        /**
         * @brief Flags that characterize the requirements for evaluating this light source.
         * @return Flags.
         */
        virtual uint32_t flags() const = 0;
    };
}
