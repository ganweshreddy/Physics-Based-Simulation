#pragma once

#include "line_fwd.hpp"

#include "base_line.hpp"

#include <vislab/core/array.hpp>
#include <vislab/core/traits.hpp>

namespace vislab
{
    /**
     * @brief Line geometry with specific internal type.
     * @tparam TArrayType Internal array type used to store vertices.
     */
    template <typename TArrayType>
    class Line : public BaseLine<TArrayType::Dimensions>
    {
        VISLAB_OBJECT(Line, BaseLine<TArrayType::Dimensions>)

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
        Line()
        {
            mVertices = std::make_shared<TArrayType>();
        }

        /**
         * @brief Copy-constructor.
         * @param other Geometry to copy from.
         */
        Line(const Line& other)
            : BaseLine<TArrayType::Dimensions>(other)
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

        /**
         * @brief Computes the length of the line.
         * @return Arc length of the line.
         */
        double length() const override
        {
            Eigen::Index numPnts = this->mVertices->getSize();
            if (numPnts <= 1)
                return typename TArrayType::Scalar(0);

            double result = 0;
            for (int i = 0; i < numPnts - 1; ++i)
            {
                const typename TArrayType::Element& a = this->mVertices->getValue(i);
                const typename TArrayType::Element& b = this->mVertices->getValue(i + 1);
                result += (a - b).norm();
            }
            return result;
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
