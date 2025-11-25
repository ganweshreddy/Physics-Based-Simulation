#pragma once

#include "surface_fwd.hpp"

#include "base_surface.hpp"

#include <vislab/core/array.hpp>
#include <vislab/core/traits.hpp>

namespace vislab
{
    /**
     * @brief Surface geometry with specific internal type.
     * @tparam TArrayType Internal array type used to store vertices.
     */
    template <typename TArrayType>
    class Surface : public BaseSurface<TArrayType::Dimensions>
    {
        VISLAB_OBJECT(Surface, BaseSurface<TArrayType::Dimensions>)

    public:
        /**
         * @brief Number of dimensions for vertex positions.
         */
        static const int Dimensions = TArrayType::Dimensions;

        /**
         * @brief Internal array type used to store the vertices.
         */
        using ArrayType = TArrayType;

        /**
         * @brief Constructor.
         */
        Surface()
        {
            mVertices = std::make_shared<TArrayType>();
        }

        /**
         * @brief Copy-constructor.
         * @param other Geometry to copy from.
         */
        Surface(const Surface& other)
            : BaseSurface<TArrayType::Dimensions>(other)
        {
            mVertices = std::shared_ptr<TArrayType>(other.mVertices->clone());
        }

        /**
         * @brief Gets the vertex position array.
         * @return Vertex position array.
         */
        std::shared_ptr<TArrayType> getVertices() { return mVertices; }

        /**
         * @brief Gets the vertex position array.
         * @return Vertex position array.
         */
        std::shared_ptr<const TArrayType> getVertices() const { return mVertices; }

        /**
         * @brief Recomputes the bounding box from the vertex buffer.
         */
        void recomputeBoundingBox() override
        {
            this->mBoundingBox.setEmpty();
            if (mVertices->getSize() == 0)
                return;
            this->mBoundingBox.extend(mVertices->getMin().template cast<double>());
            this->mBoundingBox.extend(mVertices->getMax().template cast<double>());
        }

    private:
        /**
         * @brief Gets the vertex position array.
         * @return Vertex position array.
         */
        std::shared_ptr<IArray> getVerticesImpl() override { return mVertices; }

        /**
         * @brief Gets the vertex position array.
         * @return Vertex position array.
         */
        std::shared_ptr<const IArray> getVerticesImpl() const override { return mVertices; }

        /**
         * @brief Array that stores the vertex positions.
         */
        std::shared_ptr<TArrayType> mVertices;
    };
}
