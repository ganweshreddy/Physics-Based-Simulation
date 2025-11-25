#pragma once

#include "ivector_field.hpp"

namespace vislab
{
    /**
     * @brief Base class for 3D vector fields.
     * @tparam DomainDimensions Total number of dimensions. (3=steady, 4=unsteady)
     */
    template <int64_t TDomainDimensions>
    class IVectorField3d : public IVectorField<Eigen::Vector3d, 3, TDomainDimensions>
    {
        VISLAB_INTERFACE(IVectorField3d, IVectorField<Eigen::Vector3d, 3, TDomainDimensions>)

    public:
        /**
         * @brief Type of the bounding box.
         */
        using BoundingBox = typename IVectorField<Eigen::Vector3d, 3, TDomainDimensions>::BoundingBox;

        /**
         * @brief Type of the domain coordinates.
         */
        using DomainCoord = typename IVectorField<Eigen::Vector3d, 3, TDomainDimensions>::DomainCoord;

        /**
         * @brief Type to specify the derivative degree. This type is used in the samplePartial function to specify the desired derivative of each dimension.
         */
        using DerivativeDegree = typename IVectorField<Eigen::Vector3d, 3, TDomainDimensions>::DerivativeDegree;

        /**
         * @brief Type of values in the field.
         */
        using Value = typename IVectorField<Eigen::Vector3d, 3, TDomainDimensions>::Value;

        /**
         * @brief Type of the vector-valued return types.
         */
        using Vector = typename IVectorField<Eigen::Vector3d, 3, TDomainDimensions>::Vector;

        /**
         * @brief Type of the vector-valued return types.
         */
        using Matrix = typename IVectorField<Eigen::Vector3d, 3, TDomainDimensions>::Matrix;

        /**
         * @brief Constructor.
         * @param domain Bounding box of the domain.
         */
        IVectorField3d(const BoundingBox& domain)
            : IVectorField<Eigen::Vector3d, 3, TDomainDimensions>(domain)
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
            Vector dz = sample_dxz(coord);
            Matrix result;
            result << dx, dy, dz;
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
            Vector dz = sample_dyz(coord);
            Matrix result;
            result << dx, dy, dz;
            return result;
        }

        /**
         * @brief Computes the z-partial of the Jacobian matrix.
         * @param coord Domain location at which the derivative is evaluated.
         * @return z-partial of the Jacobian matrix.
         */
        Matrix sampleJacobian_dz(const DomainCoord& coord) const
        {
            Vector dx = sample_dxz(coord);
            Vector dy = sample_dyz(coord);
            Vector dz = sample_dzz(coord);
            Matrix result;
            result << dx, dy, dz;
            return result;
        }

        /**
         * @brief Samples the x-partial derivative.
         * @param coord Domain location at which the derivative is evaluated.
         * @return x-partial derivative.
         */
        Value sample_dx(const DomainCoord& coord) const
        {
            if constexpr (TDomainDimensions == 3)
                return this->samplePartial(coord, DerivativeDegree(1, 0, 0));
            else
                return this->samplePartial(coord, DerivativeDegree(1, 0, 0, 0));
        }

        /**
         * @brief Samples the y-partial derivative.
         * @param coord Domain location at which the derivative is evaluated.
         * @return y-partial derivative.
         */
        Value sample_dy(const DomainCoord& coord) const
        {
            if constexpr (TDomainDimensions == 3)
                return this->samplePartial(coord, DerivativeDegree(0, 1, 0));
            else
                return this->samplePartial(coord, DerivativeDegree(0, 1, 0, 0));
        }

        /**
         * @brief Samples the z-partial derivative.
         * @param coord Domain location at which the derivative is evaluated.
         * @return z-partial derivative.
         */
        Value sample_dz(const DomainCoord& coord) const
        {
            if constexpr (TDomainDimensions == 3)
                return this->samplePartial(coord, DerivativeDegree(0, 0, 1));
            else
                return this->samplePartial(coord, DerivativeDegree(0, 0, 1, 0));
        }

        /**
         * @brief Samples the second x-partial derivative.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Second x-partial derivative.
         */
        Value sample_dxx(const DomainCoord& coord) const
        {
            if constexpr (TDomainDimensions == 3)
                return this->samplePartial(coord, DerivativeDegree(2, 0, 0));
            else
                return this->samplePartial(coord, DerivativeDegree(2, 0, 0, 0));
        }

        /**
         * @brief Samples the mixed xy-partial derivative.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Second mixed xy-partial derivative.
         */
        Value sample_dxy(const DomainCoord& coord) const
        {
            if constexpr (TDomainDimensions == 3)
                return this->samplePartial(coord, DerivativeDegree(1, 1, 0));
            else
                return this->samplePartial(coord, DerivativeDegree(1, 1, 0, 0));
        }

        /**
         * @brief Samples the mixed xz-partial derivative.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Second mixed xz-partial derivative.
         */
        Value sample_dxz(const DomainCoord& coord) const
        {
            if constexpr (TDomainDimensions == 3)
                return this->samplePartial(coord, DerivativeDegree(1, 0, 1));
            else
                return this->samplePartial(coord, DerivativeDegree(1, 0, 1, 0));
        }

        /**
         * @brief Samples the second y-partial derivative.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Second y-partial derivative.
         */
        Value sample_dyy(const DomainCoord& coord) const
        {
            if constexpr (TDomainDimensions == 3)
                return this->samplePartial(coord, DerivativeDegree(0, 2, 0));
            else
                return this->samplePartial(coord, DerivativeDegree(0, 2, 0, 0));
        }

        /**
         * @brief Samples the mixed yz-partial derivative.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Second mixed yz-partial derivative.
         */
        Value sample_dyz(const DomainCoord& coord) const
        {
            if constexpr (TDomainDimensions == 3)
                return this->samplePartial(coord, DerivativeDegree(0, 1, 1));
            else
                return this->samplePartial(coord, DerivativeDegree(0, 1, 1, 0));
        }

        /**
         * @brief Samples the second z-partial derivative.
         * @param coord Domain location at which the derivative is evaluated.
         * @return Second z-partial derivative.
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

        /**
         * @brief Computes the curl vector.
         * @param coord Domain location at which the quantity is evaluated.
         * @return
         */
        Vector sampleCurl(const DomainCoord& coord) const
        {
            Vector vx = sample_dx(coord);
            Vector vy = sample_dy(coord);
            Vector vz = sample_dz(coord);
            return Vector(vy.z() - vz.y(), vz.x() - vx.z(), vx.y() - vy.x());
        }

        /**
         * @brief Computes the vorticity, which is the magnitude of the curl vector.
         * @param coord Domain location at which the quantity is evaluated.
         * @return Vorticity, which is the magnitude of the curl vector.
         */
        double sampleVorticity(const DomainCoord& coord) const
        {
            return sampleCurl(coord).norm();
        }

        /**
         * @brief Computes the saddle strength.
         * @param coord Domain location at which the quantity is evaluated.
         * @return Saddle strength.
         */
        double sampleSaddleStrength(const DomainCoord& coord) const
        {
            Matrix J  = this->sampleJacobian(coord);
            auto eigs = J.eigenvalues();
            if (std::max(std::max(std::abs(eigs.imag().x()), std::abs(eigs.imag().y())), std::abs(eigs.imag().z())) >= 1E-6)
                return 0;
            Vector real(eigs.real().x(), eigs.real().y(), eigs.real().z());
            real.sortAscend();
            return -real.x() * real.z();
        }

        /**
         * @brief Computes the lambda2 value.
         * @param coord Domain location at which the quantity is evaluated.
         * @return lambda2 value.
         */
        double sampleLambda2(const DomainCoord& coord) const
        {
            Matrix J         = this->sampleJacobian(coord);
            Matrix S         = J.symmetricPart();
            Matrix O         = J.antiSymmetricPart();
            Vector eigs_real = (S * S + O * O).eigenvalues().real();
            eigs_real.sortAscend();
            return eigs_real.y();
        }

        /**
         * @brief Computes the Q vortex measure.
         * @param coord Domain location at which the quantity is evaluated.
         * @return Q vortex measure.
         */
        double sampleQHunt(const DomainCoord& coord) const
        {
            Matrix J = this->sampleJacobian(coord);
            Matrix S = J.symmetricPart();
            Matrix O = J.antiSymmetricPart();
            return (O.squaredNorm() - S.squaredNorm()) / 2;
        }

        /**
         * @brief Computes the Delta vortex criterion.
         * @param coord Domain location at which the quantity is evaluated.
         * @return Delta vortex criterion.
         */
        double sampleDelta(const DomainCoord& coord) const
        {
            Matrix J = this->sampleJacobian(coord);
            Matrix S = J.symmetricPart();
            Matrix O = J.antiSymmetricPart();
            return std::pow((O.squaredNorm() - S.squaredNorm()) / 6, 3) + std::pow(J.determinant() / 2, 2);
        }

        /**
         * @brief Computes helicity vortex measure.
         * @param coord Domain location at which the quantity is evaluated.
         * @return Helicity vortex measure.
         */
        double sampleHelicity(const DomainCoord& coord) const
        {
            Vector curl = sampleCurl(coord);
            Vector v    = this->sample(coord);
            return curl.dot(v);
        }

        /**
         * @brief Computes the gamma2 vortex measure.
         * @param coord Domain location at which the quantity is evaluated.
         * @return gamma2 vortex measure.
         */
        double sampleGamma2(const DomainCoord& coord) const
        {
            Matrix J     = this->sampleJacobian(coord);
            Matrix S     = J.symmetricPart();
            Matrix O     = J.antiSymmetricPart();
            Vector Sreal = S.eigenvalues().real(), Oreal = O.eigenvalues().real();
            return std::max(std::max(Oreal.x(), Oreal.y()), Oreal.z()) / std::max(std::max(Sreal.x(), Sreal.y()), Sreal.z());
        }
    };
}
