#pragma once

#include <vislab/field/iunsteady_vector_field2d.hpp>

#include <vislab/core/numeric_parameter.hpp>

namespace vislab
{
    /**
     * @brief Shadden's classic double gyre vector field with default parameters A=1/10, Eps=1/4, Omega=pi/5.
     */
    class DoubleGyreUVF2d : public IUnsteadyVectorField2d
    {
        VISLAB_OBJECT(DoubleGyreUVF2d, IUnsteadyVectorField2d)

    public:
        /**
         * @brief Constructor.
         * @param A Scaling factor. (default A = 0.1)
         * @param eps Unsteadiness amplitude. (default eps = 0.25)
         * @param omega Unsteadiness phase. (omega = pi / 5)
         */
        DoubleGyreUVF2d();

        /**
         * @brief Constructor.
         * @param A Scaling factor.
         * @param eps Unsteadiness amplitude.
         * @param omega Unsteadiness phase.
         */
        DoubleGyreUVF2d(double A, double eps, double omega);

        /**
         * @brief Samples the field.
         * @param coord Domain location to sample the field at.
         * @return Value at the domain location.
         */
        Eigen::Vector2d sample(const Eigen::Vector3d& coord) const override;

        /**
         * @brief Samples a partial derivative.
         * @param coord Domain location to sample the field derivative at.
         * @param dimension Enter as dimension Eigen::Vector1i(0) for the first-order x-partial, Eigen::Vector1i(1) for the first-order y-partial, Eigen::Vector2i(0,0) for the second-order x-partial, Eigen::Vector2i(0,1) for second-order mixed x- and y-partial and so on.
         * @return Field derivative at domain location.
         */
        Eigen::Vector2d samplePartial(const Eigen::Vector3d& coord, const DerivativeDegree& derivative) const override;

        /**
         * @brief Scaling factor.
         */
        DoubleParameter A;

        /**
         * @brief Unsteadiness amplitude.
         */
        DoubleParameter eps;

        /**
         * @brief Unsteadiness phase.
         */
        DoubleParameter omega;
    };
}
