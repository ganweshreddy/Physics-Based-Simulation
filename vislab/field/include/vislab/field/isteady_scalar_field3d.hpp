#pragma once

#include "iscalar_field3d.hpp"

namespace vislab
{
    /**
     * @brief Base class for 3D steady scalar fields.
     */
    class ISteadyScalarField3d : public IScalarField3d<3>
    {
        VISLAB_INTERFACE(ISteadyScalarField3d, IScalarField3d<3>)

    public:
        /**
         * @brief Type of the bounding box.
         */
        using BoundingBox = typename IScalarField3d<3>::BoundingBox;

        /**
         * @brief Type of the domain coordinates.
         */
        using DomainCoord = typename IScalarField3d<3>::DomainCoord;

        /**
         * @brief Type to specify the derivative degree. This type is used in the samplePartial function to specify the desired derivative of each dimension.
         */
        using DerivativeDegree = typename IScalarField3d<3>::DerivativeDegree;

        /**
         * @brief Type of values in the field.
         */
        using Value = typename IScalarField3d<3>::Value;

        /**
         * @brief Type of the vector-valued return types.
         */
        using Vector = typename IScalarField3d<3>::Vector;

        /**
         * @brief Type of the vector-valued return types.
         */
        using Matrix = typename IScalarField3d<3>::Matrix;

        /**
         * @brief Constructor.
         * @param domain Bounding box of the domain.
         */
        ISteadyScalarField3d(const BoundingBox& domain);
    };
}
