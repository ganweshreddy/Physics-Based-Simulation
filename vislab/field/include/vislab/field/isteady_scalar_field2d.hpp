#pragma once

#include "iscalar_field2d.hpp"

namespace vislab
{
    /**
     * @brief Base class for 2D steady scalar fields.
     */
    class ISteadyScalarField2d : public IScalarField2d<2>
    {
        VISLAB_INTERFACE(ISteadyScalarField2d, IScalarField2d<2>)

    public:
        /**
         * @brief Type of the bounding box.
         */
        using BoundingBox = typename IScalarField2d<2>::BoundingBox;

        /**
         * @brief Type of the domain coordinates.
         */
        using DomainCoord = typename IScalarField2d<2>::DomainCoord;

        /**
         * @brief Type to specify the derivative degree. This type is used in the samplePartial function to specify the desired derivative of each dimension.
         */
        using DerivativeDegree = typename IScalarField2d<2>::DerivativeDegree;

        /**
         * @brief Type of values in the field.
         */
        using Value = typename IScalarField2d<2>::Value;

        /**
         * @brief Type of the vector-valued return types.
         */
        using Vector = typename IScalarField2d<2>::Vector;

        /**
         * @brief Type of the vector-valued return types.
         */
        using Matrix = typename IScalarField2d<2>::Matrix;

        /**
         * @brief Constructor.
         * @param domain Bounding box of the domain.
         */
        ISteadyScalarField2d(const BoundingBox& domain);
    };
}
