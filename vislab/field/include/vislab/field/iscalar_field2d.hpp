#pragma once

#include "iscalar_field.hpp"

namespace vislab
{
    /**
     * @brief Base class for 2D scalar fields.
     * @tparam TDomainDimensions Total number of dimensions. (2=steady, 3=unsteady)
     */
    template <int64_t TDomainDimensions>
    class IScalarField2d : public IScalarField<2, TDomainDimensions>
    {
        VISLAB_INTERFACE(IScalarField2d, IScalarField<2, TDomainDimensions>)

    public:
        /**
         * @brief Type of the bounding box.
         */
        using BoundingBox = typename IScalarField<2, TDomainDimensions>::BoundingBox;

        /**
         * @brief Type of the domain coordinates.
         */
        using DomainCoord = typename IScalarField<2, TDomainDimensions>::DomainCoord;

        /**
         * @brief Type to specify the derivative degree. This type is used in the samplePartial function to specify the desired derivative of each dimension.
         */
        using DerivativeDegree = typename IScalarField<2, TDomainDimensions>::DerivativeDegree;

        /**
         * @brief Type of values in the field.
         */
        using Value = typename IScalarField<2, TDomainDimensions>::Value;

        /**
         * @brief Type of the vector-valued return types.
         */
        using Vector = typename IScalarField<2, TDomainDimensions>::Vector;

        /**
         * @brief Type of the vector-valued return types.
         */
        using Matrix = typename IScalarField<2, TDomainDimensions>::Matrix;

        /**
         * @brief Constructor.
         * @param domain Bounding box of the domain.
         */
        IScalarField2d(const BoundingBox& domain)
            : IScalarField<2, TDomainDimensions>(domain)
        {
        }

        /**
         * @brief Computes the co-gradient vector.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Spatial co-gradient vector.
         */
        Vector sampleCoGradient(const DomainCoord& coord) const
        {
            Vector grad = sampleGradient(coord);
            return Vector(-grad.y(), grad.x());
        }

        /**
         * @brief Samples the x-partial derivative of the field.
         * @param coord Domain location at which the derivative is evaluated.
         * @return x-partial derivative of the field.
         */
        Value sample_dx(const DomainCoord& coord) const
        {
            if constexpr (TDomainDimensions == 2)
                return this->samplePartial(coord, DerivativeDegree(1, 0));
            else
                return this->samplePartial(coord, DerivativeDegree(1, 0, 0));
        }

        /**
         * @brief Samples the y-partial derivative of the field.
         * @param coord Domain location at which the derivative is evaluated.
         * @return y-partial derivative of the field.
         */
        Value sample_dy(const DomainCoord& coord) const
        {
            if constexpr (TDomainDimensions == 2)
                return this->samplePartial(coord, DerivativeDegree(0, 1));
            else
                return this->samplePartial(coord, DerivativeDegree(0, 1, 0));
        }

        /**
         * @brief Samples the second x-partial derivative of the field.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Second x-partial derivative of the field.
         */
        Value sample_dxx(const DomainCoord& coord) const
        {
            if constexpr (TDomainDimensions == 2)
                return this->samplePartial(coord, DerivativeDegree(2, 0));
            else
                return this->samplePartial(coord, DerivativeDegree(2, 0, 0));
        }

        /**
         * @brief Samples the mixed xy-partial derivative of the field.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Second mixed xy-partial derivative of the field.
         */
        Value sample_dxy(const DomainCoord& coord) const
        {
            if constexpr (TDomainDimensions == 2)
                return this->samplePartial(coord, DerivativeDegree(1, 1));
            else
                return this->samplePartial(coord, DerivativeDegree(1, 1, 0));
        }

        /**
         * @brief Samples the second y-partial derivative of the field.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Second y-partial derivative of the field.
         */
        Value sample_dyy(const DomainCoord& coord) const
        {
            if constexpr (TDomainDimensions == 2)
                return this->samplePartial(coord, DerivativeDegree(0, 2));
            else
                return this->samplePartial(coord, DerivativeDegree(0, 2, 0));
        }

        /**
         * @brief Samples the third-order x-partial derivative.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Third-order x-partial derivative.
         */
        Value sample_dxxx(const DomainCoord& coord) const
        {
            if constexpr (TDomainDimensions == 2)
                return this->samplePartial(coord, DerivativeDegree(3, 0));
            else
                return this->samplePartial(coord, DerivativeDegree(3, 0, 0));
        }

        /**
         * @brief Samples the mixed xxy-partial derivative.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Third mixed xxy-partial derivative.
         */
        Value sample_dxxy(const DomainCoord& coord) const
        {
            if constexpr (TDomainDimensions == 2)
                return this->samplePartial(coord, DerivativeDegree(2, 1));
            else
                return this->samplePartial(coord, DerivativeDegree(2, 1, 0));
        }

        /**
         * @brief Samples the mixed xyy-partial derivative.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Third mixed xyy-partial derivative.
         */
        Value sample_dxyy(const DomainCoord& coord) const
        {
            if constexpr (TDomainDimensions == 2)
                return this->samplePartial(coord, DerivativeDegree(1, 2));
            else
                return this->samplePartial(coord, DerivativeDegree(1, 2, 0));
        }

        /**
         * @brief Samples the third-order y-partial derivative.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Third-order y-partial derivative.
         */
        Value sample_dyyy(const DomainCoord& coord) const
        {
            if constexpr (TDomainDimensions == 2)
                return this->samplePartial(coord, DerivativeDegree(0, 3));
            else
                return this->samplePartial(coord, DerivativeDegree(0, 3, 0));
        }
    };
}
