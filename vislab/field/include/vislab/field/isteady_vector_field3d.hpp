#pragma once

#include "ivector_field3d.hpp"

namespace vislab
{
    /**
     * @brief Base class for 3D steady vector fields.
     */
    class ISteadyVectorField3d : public IVectorField3d<3>
    {
        VISLAB_INTERFACE(ISteadyVectorField3d, IVectorField3d<3>)

    public:
        /**
         * @brief Type of the bounding box.
         */
        using BoundingBox = typename IVectorField3d<3>::BoundingBox;

        /**
         * @brief Type of the domain coordinates.
         */
        using DomainCoord = typename IVectorField3d<3>::DomainCoord;

        /**
         * @brief Type to specify the derivative degree. This type is used in the samplePartial function to specify the desired derivative of each dimension.
         */
        using DerivativeDegree = typename IVectorField3d<3>::DerivativeDegree;

        /**
         * @brief Type of values in the field.
         */
        using Value = typename IVectorField3d<3>::Value;

        /**
         * @brief Type of the vector-valued return types.
         */
        using Vector = typename IVectorField3d<3>::Vector;

        /**
         * @brief Type of the vector-valued return types.
         */
        using Matrix = typename IVectorField3d<3>::Matrix;

        /**
         * @brief Constructor.
         * @param domain Bounding box of the domain.
         */
        ISteadyVectorField3d(const BoundingBox& domain);

        /**
         * @brief Samples the rate of acceleration Da/Dt, also known as jerk vector.
         * @param coord Domain location at which the quantity is evaluated.
         * @return Rate of acceleration Da/Dt, also known as jerk vector.
         */
        Value sampleRateOfAcceleration(const DomainCoord& coord) const;
    };
}
