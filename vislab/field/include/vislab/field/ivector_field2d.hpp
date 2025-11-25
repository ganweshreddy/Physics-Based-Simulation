#pragma once

#include "ivector_field.hpp"

namespace vislab
{
    /**
     * @brief Base class for 2D vector fields.
     * @tparam DomainDimensions Total number of dimensions. (2=steady, 3=unsteady)
     */
    template <int64_t TDomainDimensions>
    class IVectorField2d : public IVectorField<Eigen::Vector2d, 2, TDomainDimensions>
    {
        VISLAB_INTERFACE(IVectorField2d, IVectorField<Eigen::Vector2d, 2, TDomainDimensions>)

    public:
        /**
         * @brief Type of the bounding box.
         */
        using BoundingBox = typename IVectorField<Eigen::Vector2d, 2, TDomainDimensions>::BoundingBox;

        /**
         * @brief Type of the domain coordinates.
         */
        using DomainCoord = typename IVectorField<Eigen::Vector2d, 2, TDomainDimensions>::DomainCoord;

        /**
         * @brief Type to specify the derivative degree. This type is used in the samplePartial function to specify the desired derivative of each dimension.
         */
        using DerivativeDegree = typename IVectorField<Eigen::Vector2d, 2, TDomainDimensions>::DerivativeDegree;

        /**
         * @brief Type of values in the field.
         */
        using Value = typename IVectorField<Eigen::Vector2d, 2, TDomainDimensions>::Value;

        /**
         * @brief Type of the vector-valued return types.
         */
        using Vector = typename IVectorField<Eigen::Vector2d, 2, TDomainDimensions>::Vector;

        /**
         * @brief Type of the vector-valued return types.
         */
        using Matrix = typename IVectorField<Eigen::Vector2d, 2, TDomainDimensions>::Matrix;

        /**
         * @brief Constructor.
         * @param domain Bounding box of the domain.
         */
        IVectorField2d(const BoundingBox& domain)
            : IVectorField<Eigen::Vector2d, 2, TDomainDimensions>(domain)
        {
        }

        /**
         * @brief Computes the x-partial of the Jacobian matrix.
         * @param coord Domain location at which the derivative is evaluated.
         * @return x-partial of the Jacobian matrix.
         */
        Matrix sampleJacobian_dx(const DomainCoord& coord) const
        {
            Vector dx = sample_dxx(coord);
            Vector dy = sample_dxy(coord);
            Matrix result;
            result << dx, dy;
            return result;
        }

        /**
         * @brief Computes the y-partial of the Jacobian matrix.
         * @param coord Domain location at which the derivative is evaluated.
         * @return y-partial of the Jacobian matrix.
         */
        Matrix sampleJacobian_dy(const DomainCoord& coord) const
        {
            Vector dx = sample_dxy(coord);
            Vector dy = sample_dyy(coord);
            Matrix result;
            result << dx, dy;
            return result;
        }

        /**
         * @brief Samples the signed vorticity value.
         * @param coord Domain location at which the quantity is evaluated.
         * @return Signed vorticity value.
         */
        double sampleVorticity(const DomainCoord& coord) const
        {
            Matrix O = this->sampleVorticityTensor(coord);
            return 2 * O(0, 1);
        }

        /**
         * @brief Samples the gradient of the determinant of the Jacobian.
         * @param coord Domain location at which the quantity is evaluated.
         * @return Gradient of the determinant of the Jacobian.
         */
        Vector sampleGradientDeterminantJacobian(const DomainCoord& coord) const
        {
            Vector dx  = sample_dx(coord);
            Vector dy  = sample_dy(coord);
            Vector dxx = sample_dxx(coord);
            Vector dxy = sample_dxy(coord);
            Vector dyy = sample_dyy(coord);
            double u_x = dx.x(), v_x = dx.y();
            double u_y = dy.x(), v_y = dy.y();
            double u_xx = dxx.x(), v_xx = dxx.y();
            double u_xy = dxy.x(), v_xy = dxy.y();
            double u_yy = dyy.x(), v_yy = dyy.y();
            double detJx = u_xx * v_y - u_y * v_xx + u_x * v_xy * -u_xy * v_x;
            double detJy = u_x * v_yy + u_xy * v_y - u_y * v_xy - u_yy * v_x;
            return Vector(detJx, detJy);
        }

        /**
         * @brief Samples the x-partial derivative.
         * @param coord Domain location at which the derivative is evaluated.
         * @return x-partial derivative.
         */
        Value sample_dx(const DomainCoord& coord) const
        {
            if constexpr (TDomainDimensions == 2)
                return this->samplePartial(coord, DerivativeDegree(1, 0));
            else
                return this->samplePartial(coord, DerivativeDegree(1, 0, 0));
        }

        /**
         * @brief Samples the y-partial derivative.
         * @param coord Domain location at which the derivative is evaluated.
         * @return y-partial derivative.
         */
        Value sample_dy(const DomainCoord& coord) const
        {
            if constexpr (TDomainDimensions == 2)
                return this->samplePartial(coord, DerivativeDegree(0, 1));
            else
                return this->samplePartial(coord, DerivativeDegree(0, 1, 0));
        }

        /**
         * @brief Samples the second x-partial derivative.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Second x-partial derivative.
         */
        Value sample_dxx(const DomainCoord& coord) const
        {
            if constexpr (TDomainDimensions == 2)
                return this->samplePartial(coord, DerivativeDegree(2, 0));
            else
                return this->samplePartial(coord, DerivativeDegree(2, 0, 0));
        }

        /**
         * @brief Samples the mixed xy-partial derivative.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Second mixed xy-partial derivative.
         */
        Value sample_dxy(const DomainCoord& coord) const
        {
            if constexpr (TDomainDimensions == 2)
                return this->samplePartial(coord, DerivativeDegree(1, 1));
            else
                return this->samplePartial(coord, DerivativeDegree(1, 1, 0));
        }

        /**
         * @brief Samples the second y-partial derivative.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Second y-partial derivative.
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
