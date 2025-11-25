#pragma once

#include "regular_field_fwd.hpp"

#include "base_regular_field.hpp"

#include <memory>

namespace vislab
{
    /**
     * @brief Base class for fields defined on regular grids.
     * @tparam TBaseType Field type to discretize onto a regular grid.
     * @tparam TGridType Grid type to discretize on.
     */
    template <class TBaseType, typename TGridType, typename TArrayType>
    class RegularField : public BaseRegularField<TBaseType, TGridType>
    {
        VISLAB_OBJECT(RegularField, BaseRegularField<TBaseType, TGridType>)

    public:
        /**
         * @brief Number of input dimensions.
         */
        static constexpr int Dimensions = TBaseType::Dimensions;

        /**
         * @brief Number of output dimensions.
         */
        static constexpr int Components = TBaseType::Components;

        /**
         * @brief Type of values in the field.
         */
        using Value = typename TArrayType::Element;

        /**
         * @brief Type of values in the field.
         */
        using ValueDouble = typename TBaseType::Value;

        /**
         * @brief Type of the domain coordinates.
         */
        using DomainCoord = typename TBaseType::DomainCoord;

        /**
         * @brief Type to specify the derivative degree. This type is used in the samplePartial function to specify the desired derivative of each dimension.
         */
        using DerivativeDegree = typename TBaseType::DerivativeDegree;

        /**
         * @brief Type of the underlying regular grid
         */
        using GridType = TGridType;

        /**
         * @brief Type of the underlying regular grid
         */
        using ArrayType = TArrayType;

        /**
         * @brief Type of the bounding box.
         */
        using BoundingBox = Eigen::AlignedBox<double, Dimensions>;

        /**
         * @brief Default constructor.
         */
        RegularField()
        {
        }

        /**
         * @brief Copy constructor.
         * @param other RegularField to copy from.
         */
        RegularField(const RegularField& other)
            : BaseRegularField<TBaseType, TGridType>(other)
        {
            if (other.mArray != nullptr)
                mArray = std::shared_ptr<TArrayType>(other.mArray->clone());
        }

        /**
         * @brief Gets the data array.
         * @return Underlying array that contains the data.
         */
        std::shared_ptr<TArrayType> getArray() { return mArray; }

        /**
         * @brief Gets the data array.
         * @return Underlying array that contains the data.
         */
        std::shared_ptr<const TArrayType> getArray() const { return mArray; }

        /**
         * @brief Sets the array. (The number of components must match the BaseType!)
         * @param data Array to set.
         */
        void setArray(std::shared_ptr<TArrayType> data) { mArray = data; }

        /**
         * @brief Samples the field.
         * @param coord Domain location to sample the field at.
         * @return Sampled value.
         */
        ValueDouble sample(const DomainCoord& coord) const override
        {
            return sampleImpl(coord, DerivativeDegree::Zero());
        }

        /**
         * @brief Samples a partial derivative.
         * @param coord Domain location to sample the field derivative at.
         * @param derivativeDegree Specifies the desired derivative of each dimension.
         * @return Sampled partial derivative.
         */
        ValueDouble samplePartial(const DomainCoord& coord, const DerivativeDegree& derivativeDegree) const override
        {
            return sampleImpl(coord, derivativeDegree);
        }

        /**
         * @brief Gets the vertex data at a specific grid point.
         * @param gridCoord Grid coord to get the value at.
         * @return Value at the grid coord.
         */
        Value getVertexDataAt(const typename GridType::GridCoord& gridCoord) const { return mArray->getValue(this->mGrid->getLinearIndex(gridCoord)); }

        /**
         * @brief Sets the vertex data at a specific grid point.
         * @param gridCoord Grid coord to set the value at.
         * @param value Value to set.
         */
        void setVertexDataAt(const typename GridType::GridCoord& gridCoord, const Value& value) { mArray->setValue(this->mGrid->getLinearIndex(gridCoord), value); }

        /**
         * @brief Sets the vertex data at a specific grid point if it is a scalar value.
         * @tparam D SFINAE parameter.
         * @param gridCoord Grid coord to set the value at.
         * @param value Value to set.
         * @return SFINAE return type.
         */
        template <int D = Components>
        typename std::enable_if_t<(D == 1)> setVertexDataAt(const typename GridType::GridCoord& gridCoord, const typename Value::Scalar& value)
        {
            mArray->setValue(this->mGrid->getLinearIndex(gridCoord), Value(value));
        }

    private:
        /**
         * @brief Helper function that implements the sampling of the function and its derivatives.
         * @param coord Physical domain coordinate to sample at.
         * @param derivative Specifies the desired derivative of each dimension.
         * @return Function value or its derivative.
         */
        ValueDouble sampleImpl(const DomainCoord& coord, const DerivativeDegree& derivative) const
        {
            // get some basic properties about the grid
            auto minCorner = this->mGrid->getDomain().min();
            auto maxCorner = this->mGrid->getDomain().max();
            auto res       = this->mGrid->getResolution();
            auto sep       = this->mGrid->getSpacing();

            // convert requested coordinate to index coordinates
            auto index = ((coord - minCorner).cwiseQuotient(maxCorner - minCorner)).cwiseProduct(res.template cast<double>() - DomainCoord::Ones());

            // calculate the left-most and right-most position to read from (in index coordinates) for each dimension
            Eigen::Vector<int, Dimensions> left, right;
            for (int i = 0; i < Dimensions; ++i)
            {
                int numCoefficients = 2 * ((std::max(1, derivative[i]) + 1) / 2) - 1 + this->mAccuracy;

                right[i] = (int)std::round(index[i]) + numCoefficients / 2;
                switch (this->mBoundaryBehavior[i])
                {
                case EBoundaryBehavior::Periodic:
                    right[i] = right[i] % (res[i] - 1);
                    break;
                case EBoundaryBehavior::Clamp:
                    right[i] = std::min(right[i], res[i] - 1);
                    break;
                }

                left[i] = right[i] - numCoefficients + 1;
                switch (this->mBoundaryBehavior[i])
                {
                case EBoundaryBehavior::Periodic:
                    // left bound can be out-of-range, but this is intentional to keep left<right. periodic handling is done later down below.
                    break;
                case EBoundaryBehavior::Clamp:
                    if (left[i] < 0)
                    {
                        right[i] -= left[i];
                        left[i]  = 0;
                        right[i] = std::min(right[i], res[i] - 1);
                    }
                    break;
                }
            }

            // next, we compute the coefficients per dimension.
            std::vector<double> coefficients[Dimensions];
            Eigen::Index numGridPoints = 1;
            for (int i = 0; i < Dimensions; ++i)
            {
                int numCoefficients = right[i] - left[i] + 1;
                coefficients[i].resize(numCoefficients * (derivative[i] + 1));
                calculateWeights(numCoefficients, derivative[i], sep[i], coefficients[i].data(), index[i] - left[i]);
                numGridPoints *= numCoefficients;
            }

            // lastly, we loop over the required grid points and compute the weighted sum
            ValueDouble result = ValueDouble::Zero();
            for (Eigen::Index localLinearIndex = 0; localLinearIndex < numGridPoints; ++localLinearIndex)
            {
                // compute the grid index in the stencil window
                typename GridType::GridCoord localGridIndex;
                Eigen::Index stride = 1;
                for (int d = 0; d < Dimensions - 1; ++d)
                    stride *= right[d] - left[d] + 1;
                Eigen::Index t = localLinearIndex;
                for (int64_t d = Dimensions - 1; d >= 0; --d)
                {
                    localGridIndex[d] = (int)(t / stride);
                    t                 = t % stride;
                    if (d > 0)
                        stride /= right[d - 1] - left[d - 1] + 1;
                }

                // compute the grid index in the full grid
                typename GridType::GridCoord gridIndex = localGridIndex + left;

                // do periodic handling for left bound
                for (int i = 0; i < Dimensions; ++i)
                {
                    switch (this->mBoundaryBehavior[i])
                    {
                    case EBoundaryBehavior::Periodic:
                        gridIndex[i] = gridIndex[i] % (res[i] - 1);
                        break;
                    case EBoundaryBehavior::Clamp:
                        break;
                    }
                }

                // compute the product of the weights
                double weight = 1;
                for (int i = 0; i < Dimensions; ++i)
                {
                    // we are offseting the coefficients, since the coefficients store the weights for the lower derivatives, as well.
                    weight *= coefficients[i][(right[i] - left[i] + 1) * derivative[i] + localGridIndex[i]];
                }

                // compute contribution from this stencil pixel and add to sum.
                Eigen::Index linearIndex = this->mGrid->getLinearIndex(gridIndex);
                result += mArray->getValue(linearIndex).template cast<double>() * weight;
            }

            // return the result
            return result;
        }

        /**
         * @brief Calculates weights for finite difference estimation.
         * The implementation is an adaptation of https://github.com/bjodah/finitediff (BSD license) for equispaced grids. See the paper "Generation of Finite Difference Formulas on Arbitrarily Spaced Grids", Bengt Fornberg, Mathematics of compuation, 51, 184, 1988, 699-706.
         * @tparam Real_t Scalar type for the weights.
         * @param len_g Number of sample points.
         * @param max_deriv Largest derivative to compute.
         * @param sep Grid spacing.
         * @param weights Output weights that are calculated by the function.
         * @param around Location in grid coordinates (zero-based) at which the derivative is evaluated.
         */
        template <typename Real_t>
        void calculateWeights(const unsigned len_g, const unsigned max_deriv, const double sep, Real_t* const weights, const Real_t around = 0) const
        {
            // if not enough sample points, return zero.
            if (len_g < max_deriv + 1)
            {
                for (unsigned i = 0; i < len_g * (max_deriv + 1); ++i)
                    weights[i] = 0; // clear weights
                return;
            }
            Real_t c1, c4, c5;
            c1         = 1;
            c4         = (0 - around) * sep;
            weights[0] = 1;
            for (unsigned i = 1; i < len_g * (max_deriv + 1); ++i)
                weights[i] = 0; // clear weights
            for (unsigned i = 1; i < len_g; ++i)
            {
                const int mn = std::min(i, max_deriv);
                Real_t c2    = 1;
                c5           = c4;
                c4           = (i - around) * sep;
                for (unsigned j = 0; j < i; ++j)
                {
                    const Real_t c3   = (i - j) * sep;
                    const Real_t c3_r = 1 / c3;
                    c2                = c2 * c3;
                    if (j == i - 1)
                    {
                        const Real_t c2_r = 1 / c2;
                        for (int k = mn; k >= 1; --k)
                        {
                            const Real_t tmp1      = weights[i - 1 + (k - 1) * len_g];
                            const Real_t tmp2      = weights[i - 1 + k * len_g];
                            weights[i + k * len_g] = c1 * (k * tmp1 - c5 * tmp2) * c2_r;
                        }
                        weights[i] = -c1 * c5 * weights[i - 1] * c2_r;
                    }
                    for (unsigned k = mn; k >= 1; --k)
                    {
                        const Real_t tmp1      = weights[j + k * len_g];
                        const Real_t tmp2      = weights[j + (k - 1) * len_g];
                        weights[j + k * len_g] = (c4 * tmp1 - k * tmp2) * c3_r;
                    }
                    weights[j] = c4 * weights[j] * c3_r;
                }
                c1 = c2;
            }
        }

        /**
         * @brief Gets the data array.
         * @return Underlying array that contains the data.
         */
        virtual std::shared_ptr<typename BaseRegularField<TBaseType, TGridType>::ArrayType> getArrayImpl() { return mArray; }

        /**
         * @brief Gets the data array.
         * @return Underlying array that contains the data.
         */
        virtual std::shared_ptr<const typename BaseRegularField<TBaseType, TGridType>::ArrayType> getArrayImpl() const { return mArray; }

        /**
         * @brief Gets the vertex data at a specific grid point.
         * @param gridCoord Grid coord to get the value at.
         * @return Value at the grid coord.
         */
        virtual ValueDouble getVertexDataAtDoubleImpl(const typename GridType::GridCoord& gridCoord) const
        {
            return mArray->getValueDouble(this->mGrid->getLinearIndex(gridCoord));
        }

        /**
         * @brief Sets the vertex data at a specific grid point.
         * @param gridCoord Grid coord to set the value at.
         * @param value Value to set.
         */
        virtual void setVertexDataAtDoubleImpl(const typename GridType::GridCoord& gridCoord, const ValueDouble& value)
        {
            return mArray->setValueDouble(this->mGrid->getLinearIndex(gridCoord), value);
        }

        /**
         * @brief Data that is stored on the grid points.
         */
        std::shared_ptr<TArrayType> mArray;
    };
}
