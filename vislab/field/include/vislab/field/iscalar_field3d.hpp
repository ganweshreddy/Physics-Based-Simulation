#pragma once

#include "iscalar_field.hpp"

namespace vislab
{
    /**
     * @brief Base class for 3D scalar fields.
     * @tparam TDomainDimensions Total number of dimensions. (3=steady, 4=unsteady)
     */
    template <int64_t TDomainDimensions>
    class IScalarField3d : public IScalarField<3, TDomainDimensions>
    {
        VISLAB_INTERFACE(IScalarField3d, IScalarField<3, TDomainDimensions>)

    public:
        /**
         * @brief Type of the bounding box.
         */
        using BoundingBox = typename IScalarField<3, TDomainDimensions>::BoundingBox;

        /**
         * @brief Type of the domain coordinates.
         */
        using DomainCoord = typename IScalarField<3, TDomainDimensions>::DomainCoord;

        /**
         * @brief Type to specify the derivative degree. This type is used in the samplePartial function to specify the desired derivative of each dimension.
         */
        using DerivativeDegree = typename IScalarField<3, TDomainDimensions>::DerivativeDegree;

        /**
         * @brief Type of values in the field.
         */
        using Value = typename IScalarField<3, TDomainDimensions>::Value;

        /**
         * @brief Type of the vector-valued return types.
         */
        using Vector = typename IScalarField<3, TDomainDimensions>::Vector;

        /**
         * @brief Type of the vector-valued return types.
         */
        using Matrix = typename IScalarField<3, TDomainDimensions>::Matrix;

        /**
         * @brief Constructor.
         * @param domain Bounding box of the domain.
         */
        IScalarField3d(const BoundingBox& domain)
            : IScalarField<3, TDomainDimensions>(domain)
        {
        }

        /**
         * @brief Samples the x-partial derivative of the field.
         * @param coord Domain location at which the derivative is evaluated.
         * @return x-partial derivative of the field.
         */
        Value sample_dx(const DomainCoord& coord) const
        {
            if constexpr (TDomainDimensions == 3)
                return this->samplePartial(coord, DerivativeDegree(1, 0, 0));
            else
                return this->samplePartial(coord, DerivativeDegree(1, 0, 0, 0));
        }

        /**
         * @brief Samples the y-partial derivative of the field.
         * @param coord Domain location at which the derivative is evaluated.
         * @return y-partial derivative of the field.
         */
        Value sample_dy(const DomainCoord& coord) const
        {
            if constexpr (TDomainDimensions == 3)
                return this->samplePartial(coord, DerivativeDegree(0, 1, 0));
            else
                return this->samplePartial(coord, DerivativeDegree(0, 1, 0, 0));
        }

        /**
         * @brief Samples the z-partial derivative of the field.
         * @param coord Domain location at which the derivative is evaluated.
         * @return z-partial derivative of the field.
         */
        Value sample_dz(const DomainCoord& coord) const
        {
            if constexpr (TDomainDimensions == 3)
                return this->samplePartial(coord, DerivativeDegree(0, 0, 1));
            else
                return this->samplePartial(coord, DerivativeDegree(0, 0, 1, 0));
        }

        /**
         * @brief Samples the second x-partial derivative of the field.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Second x-partial derivative of the field.
         */
        Value sample_dxx(const DomainCoord& coord) const
        {
            if constexpr (TDomainDimensions == 3)
                return this->samplePartial(coord, DerivativeDegree(2, 0, 0));
            else
                return this->samplePartial(coord, DerivativeDegree(2, 0, 0, 0));
        }

        /**
         * @brief Samples the mixed xy-partial derivative of the field.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Second mixed xy-partial derivative of the field.
         */
        Value sample_dxy(const DomainCoord& coord) const
        {
            if constexpr (TDomainDimensions == 3)
                return this->samplePartial(coord, DerivativeDegree(1, 1, 0));
            else
                return this->samplePartial(coord, DerivativeDegree(1, 1, 0, 0));
        }

        /**
         * @brief Samples the mixed xz-partial derivative of the field.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Second mixed xz-partial derivative of the field.
         */
        Value sample_dxz(const DomainCoord& coord) const
        {
            if constexpr (TDomainDimensions == 3)
                return this->samplePartial(coord, DerivativeDegree(1, 0, 1));
            else
                return this->samplePartial(coord, DerivativeDegree(1, 0, 1, 0));
        }

        /**
         * @brief Samples the second y-partial derivative of the field.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Second y-partial derivative of the field.
         */
        Value sample_dyy(const DomainCoord& coord) const
        {
            if constexpr (TDomainDimensions == 3)
                return this->samplePartial(coord, DerivativeDegree(0, 2, 0));
            else
                return this->samplePartial(coord, DerivativeDegree(0, 2, 0, 0));
        }

        /**
         * @brief Samples the mixed yz-partial derivative of the field.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Second mixed yz-partial derivative of the field.
         */
        Value sample_dyz(const DomainCoord& coord) const
        {
            if constexpr (TDomainDimensions == 3)
                return this->samplePartial(coord, DerivativeDegree(0, 1, 1));
            else
                return this->samplePartial(coord, DerivativeDegree(0, 1, 1, 0));
        }

        /**
         * @brief Samples the second z-partial derivative of the field.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Second z-partial derivative of the field.
         */
        Value sample_dzz(const DomainCoord& coord) const
        {
            if constexpr (TDomainDimensions == 3)
                return this->samplePartial(coord, DerivativeDegree(0, 0, 2));
            else
                return this->samplePartial(coord, DerivativeDegree(0, 0, 2, 0));
        }

        /**
         * @brief Samples the third-order x-partial derivative.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Third-order x-partial derivative.
         */
        Value sample_dxxx(const DomainCoord& coord) const
        {
            if constexpr (TDomainDimensions == 3)
                return this->samplePartial(coord, DerivativeDegree(3, 0, 0));
            else
                return this->samplePartial(coord, DerivativeDegree(3, 0, 0, 0));
        }

        /**
         * @brief Samples the mixed xxy-partial derivative.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Third mixed xxy-partial derivative.
         */
        Value sample_dxxy(const DomainCoord& coord) const
        {
            if constexpr (TDomainDimensions == 3)
                return this->samplePartial(coord, DerivativeDegree(2, 1, 0));
            else
                return this->samplePartial(coord, DerivativeDegree(2, 1, 0, 0));
        }

        /**
         * @brief Samples the mixed xxz-partial derivative.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Third mixed xxz-partial derivative.
         */
        Value sample_dxxz(const DomainCoord& coord) const
        {
            if constexpr (TDomainDimensions == 3)
                return this->samplePartial(coord, DerivativeDegree(2, 0, 1));
            else
                return this->samplePartial(coord, DerivativeDegree(2, 0, 1, 0));
        }

        /**
         * @brief Samples the mixed xyy-partial derivative.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Third mixed xyy-partial derivative.
         */
        Value sample_dxyy(const DomainCoord& coord) const
        {
            if constexpr (TDomainDimensions == 3)
                return this->samplePartial(coord, DerivativeDegree(1, 2, 0));
            else
                return this->samplePartial(coord, DerivativeDegree(1, 2, 0, 0));
        }

        /**
         * @brief Samples the mixed xyz-partial derivative.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Third mixed xyz-partial derivative.
         */
        Value sample_dxyz(const DomainCoord& coord) const
        {
            if constexpr (TDomainDimensions == 3)
                return this->samplePartial(coord, DerivativeDegree(1, 1, 1));
            else
                return this->samplePartial(coord, DerivativeDegree(1, 1, 1, 0));
        }

        /**
         * @brief Samples the mixed xzz-partial derivative.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Third mixed xzz-partial derivative.
         */
        Value sample_dxzz(const DomainCoord& coord) const
        {
            if constexpr (TDomainDimensions == 3)
                return this->samplePartial(coord, DerivativeDegree(1, 0, 2));
            else
                return this->samplePartial(coord, DerivativeDegree(1, 0, 2, 0));
        }

        /**
         * @brief Samples the third-order y-partial derivative.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Third-order y-partial derivative.
         */
        Value sample_dyyy(const DomainCoord& coord) const
        {
            if constexpr (TDomainDimensions == 3)
                return this->samplePartial(coord, DerivativeDegree(0, 3, 0));
            else
                return this->samplePartial(coord, DerivativeDegree(0, 3, 0, 0));
        }

        /**
         * @brief Samples the mixed yyz-partial derivative.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Third mixed yyz-partial derivative.
         */
        Value sample_dyyz(const DomainCoord& coord) const
        {
            if constexpr (TDomainDimensions == 3)
                return this->samplePartial(coord, DerivativeDegree(0, 2, 1));
            else
                return this->samplePartial(coord, DerivativeDegree(0, 2, 1, 0));
        }

        /**
         * @brief Samples the mixed yzz-partial derivative.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Third mixed yzz-partial derivative.
         */
        Value sample_dyzz(const DomainCoord& coord) const
        {
            if constexpr (TDomainDimensions == 3)
                return this->samplePartial(coord, DerivativeDegree(0, 1, 2));
            else
                return this->samplePartial(coord, DerivativeDegree(0, 1, 2, 0));
        }

        /**
         * @brief Samples the third-order z-partial derivative.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Third-order z-partial derivative.
         */
        Value sample_dzzz(const DomainCoord& coord) const
        {
            if constexpr (TDomainDimensions == 3)
                return this->samplePartial(coord, DerivativeDegree(0, 0, 3));
            else
                return this->samplePartial(coord, DerivativeDegree(0, 0, 3, 0));
        }
    };
}
