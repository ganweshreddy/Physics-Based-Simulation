#pragma once

#include "ivector_field2d.hpp"

namespace vislab
{
    /**
     * @brief Base class for 2D steady vector fields.
     */
    class ISteadyVectorField2d : public IVectorField2d<2>
    {
        VISLAB_INTERFACE(ISteadyVectorField2d, IVectorField2d<2>)

    public:
        /**
         * @brief Type of the bounding box.
         */
        using BoundingBox = typename IVectorField2d<2>::BoundingBox;

        /**
         * @brief Type of the domain coordinates.
         */
        using DomainCoord = typename IVectorField2d<2>::DomainCoord;

        /**
         * @brief Type to specify the derivative degree. This type is used in the samplePartial function to specify the desired derivative of each dimension.
         */
        using DerivativeDegree = typename IVectorField2d<2>::DerivativeDegree;

        /**
         * @brief Type of values in the field.
         */
        using Value = typename IVectorField2d<2>::Value;

        /**
         * @brief Type of the vector-valued return types.
         */
        using Vector = typename IVectorField2d<2>::Vector;

        /**
         * @brief Type of the vector-valued return types.
         */
        using Matrix = typename IVectorField2d<2>::Matrix;

        /**
         * @brief Constructor.
         * @param domain Bounding box of the domain.
         */
        ISteadyVectorField2d(const BoundingBox& domain);

        /**
         * @brief Samples the rate of acceleration Da/Dt, also known as jerk vector.
         * @param coord Domain location at which the quantity is evaluated.
         * @return Rate of acceleration Da/Dt, also known as jerk vector.
         */
        Value sampleRateOfAcceleration(const DomainCoord& coord) const;
    };
}
