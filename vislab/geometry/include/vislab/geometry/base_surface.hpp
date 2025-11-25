#pragma once

#include "base_surface_fwd.hpp"

#include "isurface.hpp"

#include <vislab/core/base_array.hpp>

namespace vislab
{
    /**
     * @brief Base class for surface geometry with a certain dimensionality.
     * @tparam TDimensions Number of dimensions.
     */
    template <int64_t TDimensions>
    class BaseSurface : public ISurface
    {
        VISLAB_INTERFACE(BaseSurface, ISurface)

    public:
        /**
         * @brief Number of spatial dimensions.
         */
        static const int64_t Dimensions = TDimensions;

        /**
         * @brief Bounding box type for this geometry.
         */
        using BoundingBox = Eigen::AlignedBox<double, TDimensions>;

        /**
         * @brief Constructor.
         */
        BaseSurface()
        {
            mBoundingBox.setEmpty();
        }

        /**
         * @brief Copy-constructor.
         * @param other Geometry to copy from.
         */
        BaseSurface(const BaseSurface& other)
            : ISurface(other)
        {
            mBoundingBox = other.mBoundingBox;
        }

        /**
         * @brief Gets the vertex position array.
         * @return Vertex position array.
         */
        std::shared_ptr<BaseArray<TDimensions>> getVertices() { return std::static_pointer_cast<BaseArray<TDimensions>>(getVerticesImpl()); }

        /**
         * @brief Gets the vertex position array.
         * @return Vertex position array.
         */
        std::shared_ptr<const BaseArray<TDimensions>> getVertices() const { return std::static_pointer_cast<const BaseArray<TDimensions>>(getVerticesImpl()); }

        /**
         * @brief Gets the bounding box of the vertices. Note that recomputeBoundingBox has to be called first!
         * @return Bounding box of this geometry.
         */
        const BoundingBox& getBoundingBox() const { return mBoundingBox; }

        /**
         * @brief Tests if two point geometries are equal.
         * @param other Points to compare to.
         * @return True, if the geometries are equal.
         */
        bool isEqual(const ISurface* other) const override
        {
            if (!ISurface::isEqual(other))
                return false;
            const BaseSurface* otherTyped = dynamic_cast<const BaseSurface*>(other);
            if (otherTyped == nullptr)
                return false;
            if (!mBoundingBox.isApprox(otherTyped->mBoundingBox))
                return false;
            return true;
        }

        /**
         * @brief Removes all vertices and attributes.
         */
        void clear() override
        {
            ISurface::clear();
            mBoundingBox.setEmpty();
        }

    protected:
        /**
         * @brief Bounding box of the vertices.
         */
        BoundingBox mBoundingBox;
    };
}
