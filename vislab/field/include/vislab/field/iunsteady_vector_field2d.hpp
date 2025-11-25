#pragma once

#include "ivector_field2d.hpp"

namespace vislab
{
    /**
     * @brief Base class for 2D unsteady vector fields.
     */
    class IUnsteadyVectorField2d : public IVectorField2d<3>
    {
        VISLAB_INTERFACE(IUnsteadyVectorField2d, IVectorField2d<3>)

    public:
        /**
         * @brief Type of the bounding box.
         */
        using BoundingBox = typename IVectorField2d<3>::BoundingBox;

        /**
         * @brief Type of the domain coordinates.
         */
        using DomainCoord = typename IVectorField2d<3>::DomainCoord;

        /**
         * @brief Type to specify the derivative degree. This type is used in the samplePartial function to specify the desired derivative of each dimension.
         */
        using DerivativeDegree = typename IVectorField2d<3>::DerivativeDegree;

        /**
         * @brief Type of values in the field.
         */
        using Value = typename IVectorField2d<3>::Value;

        /**
         * @brief Type of the vector-valued return types.
         */
        using Vector = typename IVectorField2d<3>::Vector;

        /**
         * @brief Type of the vector-valued return types.
         */
        using Matrix = typename IVectorField2d<3>::Matrix;

        /**
         * @brief Constructor.
         * @param domain Bounding box of the domain.
         */
        IUnsteadyVectorField2d(const BoundingBox& domain);

        /**
         * @brief Samples the t-partial derivative.
         * @param coord Domain location at which the derivative is evaluated.
         * @return t-partial derivative.
         */
        Value sample_dt(const DomainCoord& coord) const;

        /**
         * @brief Samples the mixed xt-partial derivative.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Second mixed xt-partial derivative.
         */
        Value sample_dxt(const DomainCoord& coord) const;

        /**
         * @brief Samples the mixed yt-partial derivative.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Second mixed yt-partial derivative.
         */
        Value sample_dyt(const DomainCoord& coord) const;

        /**
         * @brief Samples the mixed tt-partial derivative.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Second mixed tt-partial derivative.
         */
        Value sample_dtt(const DomainCoord& coord) const;

        /**
         * @brief Samples the mixed xxt-partial derivative.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Third mixed xxt-partial derivative.
         */
        Value sample_dxxt(const DomainCoord& coord) const;

        /**
         * @brief Samples the mixed xyt-partial derivative.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Third mixed xyt-partial derivative.
         */
        Value sample_dxyt(const DomainCoord& coord) const;

        /**
         * @brief Samples the mixed yyt-partial derivative.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Third mixed yyt-partial derivative.
         */
        Value sample_dyyt(const DomainCoord& coord) const;

        /**
         * @brief Samples the mixed xtt-partial derivative.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Third mixed xtt-partial derivative.
         */
        Value sample_dxtt(const DomainCoord& coord) const;

        /**
         * @brief Samples the mixed ytt-partial derivative.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Third mixed ytt-partial derivative.
         */
        Value sample_dytt(const DomainCoord& coord) const;

        /**
         * @brief Samples the third-order t-partial derivative.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Third-order t-partial derivative.
         */
        Value sample_dttt(const DomainCoord& coord) const;

        /**
         * @brief Computes the t-partial of the Jacobian matrix.
         * @param coord Domain location at which the derivative is evaluated.
         * @return t-partial of the Jacobian matrix.
         */
        Matrix sampleJacobian_dt(const DomainCoord& coord) const;

        /**
         * @brief Samples the material derivative of the Jacobian.
         * @param coord Domain location at which the quantity is evaluated.
         * @return Material derivative of the Jacobian.
         */
        Matrix sampleMaterialDerivativeJacobian(const DomainCoord& coord) const;

        /**
         * @brief Samples the rate of acceleration Da/Dt, also known as jerk vector.
         * @param coord Domain location at which the quantity is evaluated.
         * @return Rate of acceleration Da/Dt, also known as jerk vector.
         */
        Value sampleRateOfAcceleration(const DomainCoord& coord) const;

        /**
         * @brief Samples the feature flow field: -J^-1 * v_t.
         * @param coord Domain location at which the quantity is evaluated.
         * @return Feature flow field.
         */
        Value sampleFeatureFlowField(const DomainCoord& coord) const;
    };
}
