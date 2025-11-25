#pragma once

#include "base_field.hpp"

namespace vislab
{
    /**
     * @brief Base class for scalar fields in arbitrary dimension.
     * @tparam SpatialDimensions Number of spatial dimensions.
     * @tparam DomainDimensions Total number of dimensions.
     */
    template <int64_t SpatialDimensions, int64_t DomainDimensions>
    class IScalarField : public BaseField<Eigen::Vector1d, SpatialDimensions, DomainDimensions>
    {
        VISLAB_INTERFACE(IScalarField, BaseField<Eigen::Vector1d, SpatialDimensions, DomainDimensions>)

    public:
        /**
         * @brief Type of the bounding box.
         */
        using BoundingBox = typename BaseField<Eigen::Vector1d, SpatialDimensions, DomainDimensions>::BoundingBox;

        /**
         * @brief Type of the domain coordinates.
         */
        using DomainCoord = typename BaseField<Eigen::Vector1d, SpatialDimensions, DomainDimensions>::DomainCoord;

        /**
         * @brief Type to specify the derivative degree. This type is used in the samplePartial function to specify the desired derivative of each dimension.
         */
        using DerivativeDegree = typename BaseField<Eigen::Vector1d, SpatialDimensions, DomainDimensions>::DerivativeDegree;

        /**
         * @brief Type of values in the field.
         */
        using Value = typename BaseField<Eigen::Vector1d, SpatialDimensions, DomainDimensions>::Value;

        /**
         * @brief Type of the vector-valued return types.
         */
        using Vector = Eigen::Matrix<double, SpatialDimensions, 1>;

        /**
         * @brief Type of the vector-valued return types.
         */
        using Matrix = Eigen::Matrix<double, SpatialDimensions, SpatialDimensions>;

        /**
         * @brief Constructor.
         * @param domain Bounding box of the domain.
         */
        IScalarField(const BoundingBox& domain)
            : BaseField<Eigen::Vector1d, SpatialDimensions, DomainDimensions>(domain)
        {
        }

        /**
         * @brief Computes the gradient vector.
         * @param coord Domain location.
         * @return Spatial gradient at domain location.
         */
        Vector sampleGradient(const DomainCoord& coord) const
        {
            Vector gradient;
            for (int i = 0; i < SpatialDimensions; ++i)
            {
                DerivativeDegree derivative = DerivativeDegree::Zero();
                derivative[i]               = 1;
                gradient[i]                 = this->samplePartial(coord, derivative).x();
            }
            return gradient;
        }

        /**
         * @brief Computes the Hessian matrix.
         * @param coord Domain location.
         * @return Spatial hessian at domain location.
         */
        Matrix sampleHessian(const DomainCoord& coord) const
        {
            Matrix hessian;
            for (int r = 0; r < SpatialDimensions; ++r)
                for (int c = r; c < SpatialDimensions; ++c)
                {
                    hessian(r, c) = this->samplePartial(coord, Eigen::Vector2i(r, c)).x();
                    if (r != c)
                        hessian(c, r) = hessian(r, c);
                }
            return hessian;
        }
    };
}
