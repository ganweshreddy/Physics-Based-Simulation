#pragma once

#include "iscalar_field3d.hpp"

namespace vislab
{
    /**
     * @brief Base class for 3D unsteady scalar fields.
     */
    class IUnsteadyScalarField3d : public IScalarField3d<4>
    {
        VISLAB_INTERFACE(IUnsteadyScalarField3d, IScalarField3d<4>)

    public:
        /**
         * @brief Type of the bounding box.
         */
        using BoundingBox = typename IScalarField3d<4>::BoundingBox;

        /**
         * @brief Type of the domain coordinates.
         */
        using DomainCoord = typename IScalarField3d<4>::DomainCoord;

        /**
         * @brief Type to specify the derivative degree. This type is used in the samplePartial function to specify the desired derivative of each dimension.
         */
        using DerivativeDegree = typename IScalarField3d<4>::DerivativeDegree;

        /**
         * @brief Type of values in the field.
         */
        using Value = typename IScalarField3d<4>::Value;

        /**
         * @brief Type of the vector-valued return types.
         */
        using Vector = typename IScalarField3d<4>::Vector;

        /**
         * @brief Type of the vector-valued return types.
         */
        using Matrix = typename IScalarField3d<4>::Matrix;

        /**
         * @brief Constructor.
         * @param domain
         */
        IUnsteadyScalarField3d(const BoundingBox& domain);

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
         * @brief Samples the mixed zt-partial derivative.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Second mixed zt-partial derivative.
         */
        Value sample_dzt(const DomainCoord& coord) const;

        /**
         * @brief Samples the second t-partial derivative.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Second t-partial derivative.
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
         * @brief Samples the mixed xzt-partial derivative.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Third mixed xzt-partial derivative.
         */
        Value sample_dxzt(const DomainCoord& coord) const;

        /**
         * @brief Samples the mixed yyt-partial derivative.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Third mixed yyt-partial derivative.
         */
        Value sample_dyyt(const DomainCoord& coord) const;

        /**
         * @brief Samples the mixed yzt-partial derivative.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Third mixed yzt-partial derivative.
         */
        Value sample_dyzt(const DomainCoord& coord) const;

        /**
         * @brief Samples the mixed zzt-partial derivative.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Third mixed zzt-partial derivative.
         */
        Value sample_dzzt(const DomainCoord& coord) const;

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
         * @brief Samples the mixed ztt-partial derivative.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Third mixed ztt-partial derivative.
         */
        Value sample_dztt(const DomainCoord& coord) const;

        /**
         * @brief Samples the third-order t-partial derivative.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Third-order t-partial derivative.
         */
        Value sample_dttt(const DomainCoord& coord) const;
    };
}
