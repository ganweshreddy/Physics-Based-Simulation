#pragma once

#include "base_regular_field_fwd.hpp"

#include "isteady_scalar_field2d.hpp"
#include "isteady_scalar_field3d.hpp"
#include "isteady_vector_field2d.hpp"
#include "isteady_vector_field3d.hpp"
#include "iunsteady_scalar_field2d.hpp"
#include "iunsteady_scalar_field3d.hpp"
#include "iunsteady_vector_field2d.hpp"
#include "iunsteady_vector_field3d.hpp"
#include "regular_grid.hpp"

#include <vislab/core/array.hpp>

#include <memory>

namespace vislab
{
    /**
     * @brief Enumeration of boundary behaviors supported by the regular field.
     */
    enum class EBoundaryBehavior
    {
        /**
         * @brief Samples outside the domain get clamped to the interior.
         */
        Clamp,

        /**
         * @brief The domain is periodic and samples that move out on one side are consided to be coming in from the opposite side.
         */
        Periodic
    };

    /**
     * @brief Base class for fields defined on regular grids.
     * @tparam TBaseType Field type to discretize onto a regular grid.
     * @tparam TGridType Grid type to discretize on.
     */
    template <class TBaseType, typename TGridType>
    class BaseRegularField : public TBaseType
    {
        VISLAB_INTERFACE(BaseRegularField, TBaseType)

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
        using Value = typename TBaseType::Value;

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
        using ArrayType = BaseArray<Components>;

        /**
         * @brief Type of the bounding box.
         */
        using BoundingBox = Eigen::AlignedBox<double, Dimensions>;

        /**
         * @brief Typeof the values in double precision.
         */
        using ValueDouble = Eigen::Vector<double, ArrayType::Dimensions>;

        /**
         * @brief Default constructor.
         */
        BaseRegularField()
            : TBaseType(BoundingBox())
            , mAccuracy(2)
        {
            for (int i = 0; i < Dimensions; ++i)
                mBoundaryBehavior[i] = EBoundaryBehavior::Clamp;
        }

        /**
         * @brief Copy constructor.
         * @param other RegularField to copy from.
         */
        BaseRegularField(const BaseRegularField& other)
            : TBaseType(BoundingBox())
            , mBoundaryBehavior(other.mBoundaryBehavior)
            , mAccuracy(other.mAccuracy)
        {
            if (other.mGrid != nullptr)
                mGrid = std::shared_ptr<GridType>(other.mGrid->clone());
        }

        /**
         * @brief Sets the grid.
         * @param grid Grid to set.
         */
        void setGrid(std::shared_ptr<GridType> grid) { mGrid = grid; }

        ///**
        // * @brief Sets the array. (The number of components must match the BaseType!)
        // * @param data Array to set.
        // */
        // void setArray(std::shared_ptr<ArrayType> data) { mArray = data; }

        /**
         * @brief Gets the grid.
         * @return Underlying grid.
         */
        std::shared_ptr<const GridType> getGrid() const { return mGrid; }

        /**
         * @brief Gets the data array.
         * @return Underlying array that contains the data.
         */
        std::shared_ptr<const ArrayType> getArray() const { return getArrayImpl(); }

        /**
         * @brief Gets the grid.
         * @return Underlying grid.
         */
        std::shared_ptr<GridType> getGrid() { return mGrid; }

        /**
         * @brief Gets the data array.
         * @return Underlying array that contains the data.
         */
        std::shared_ptr<ArrayType> getArray() { return getArrayImpl(); }

        /**
         * @brief Gets the vertex data at a specific grid point.
         * @param gridCoord Grid coord to get the value at.
         * @return Value at the grid coord.
         */
        ValueDouble getVertexDataAtDouble(const typename GridType::GridCoord& gridCoord) const { return getVertexDataAtDoubleImpl(gridCoord); }

        /**
         * @brief Sets the vertex data at a specific grid point.
         * @param gridCoord Grid coord to set the value at.
         * @param value Value to set.
         */
        void setVertexDataAtDouble(const typename GridType::GridCoord& gridCoord, const ValueDouble& value) { setVertexDataAtDoubleImpl(gridCoord, value); }

        /**
         * @brief Gets the domain of the field.
         * @return Domain bounding box.
         */
        const BoundingBox& getDomain() const override { return mGrid->getDomain(); }

        /**
         * @brief Sets the domain of the field.
         * @param domain Domain bounding box to set.
         */
        void setDomain(const BoundingBox& domain) override
        {
            this->mDomain = domain;
            mGrid->setDomain(domain);
        }

        /**
         * @brief Gets the boundary behavior for a certain dimension.
         * @param dimension Dimension to get boundary behavior for.
         * @return Boundary behavior of requested dimension.
         */
        const EBoundaryBehavior& getBoundaryBehavior(const int& dimension) const
        {
            return mBoundaryBehavior[dimension];
        }

        /**
         * @brief Sets the boundary behavior for a certain dimension.
         * @param dimension Dimension to set boundary behavior for.
         * @param behavior Boundary behavior to set.
         */
        void setBoundaryBehavior(const int& dimension, const EBoundaryBehavior& behavior)
        {
            mBoundaryBehavior[dimension] = behavior;
        }

        /**
         * @brief Gets the accuracy of the Newton interpolation. This parameter determines the degree of the interpolation polynomial.
         * @return Accuracy of Newton interpolation.
         */
        const int& getAccuracy() const
        {
            return mAccuracy;
        }

        /**
         * @brief Sets the accuracy of the Newton interpolation. This parameter determines the degree of the interpolation polynomial.
         * @param accuracy Accuracy of Newton interpolation to set.
         */
        void setAccuracy(const int& accuracy)
        {
            mAccuracy = accuracy;
        }

    protected:
        /**
         * @brief Gets the data array.
         * @return Underlying array that contains the data.
         */
        virtual std::shared_ptr<ArrayType> getArrayImpl() = 0;

        /**
         * @brief Gets the data array.
         * @return Underlying array that contains the data.
         */
        virtual std::shared_ptr<const ArrayType> getArrayImpl() const = 0;

        /**
         * @brief Gets the vertex data at a specific grid point.
         * @param gridCoord Grid coord to get the value at.
         * @return Value at the grid coord.
         */
        virtual ValueDouble getVertexDataAtDoubleImpl(const typename GridType::GridCoord& gridCoord) const = 0;

        /**
         * @brief Sets the vertex data at a specific grid point.
         * @param gridCoord Grid coord to set the value at.
         * @param value Value to set.
         */
        virtual void setVertexDataAtDoubleImpl(const typename GridType::GridCoord& gridCoord, const ValueDouble& value) = 0;

        /**
         * @brief Grid on which the data is stored.
         */
        std::shared_ptr<GridType> mGrid;

        /**
         * @brief Sets the behavior on the axis-aligned boundaries for samples that land outside of the domain.
         */
        std::array<EBoundaryBehavior, Dimensions> mBoundaryBehavior;

        /**
         * @brief Accuracy of the Newton interpolation. This parameter determines the degree of the interpolation polynomial.
         */
        int mAccuracy;
    };
}
