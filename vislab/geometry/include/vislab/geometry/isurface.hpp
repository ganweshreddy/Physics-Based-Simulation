#pragma once

#include <memory>
#include <vislab/core/base_array_fwd.hpp>
#include <vislab/core/data.hpp>
#include <vislab/core/types.hpp>

namespace vislab
{
    class IArray;
    class Attributes;

    /**
     * @brief Interface for class that stores an individual surface.
     */
    class ISurface : public Data
    {
        VISLAB_INTERFACE(ISurface, Data)

    public:
        /**
         * @brief Primitive topology.
         */
        enum class ETopology
        {
            /**
             * @brief List of triangles.
             */
            TriangleList,

            /**
             * @brief Strip of triangles.
             */
            TriangleStrip,

            /**
             * @brief List of triangular patches.
             */
            PatchList3,

            /**
             * @brief List of quad patches.
             */
            PatchList4,
        };

        /**
         * @brief Constructor.
         */
        ISurface();

        /**
         * @brief Copy-constructor.
         * @param other Surface to copy from.
         */
        ISurface(const ISurface& other);

        /**
         * @brief Destructor.
         */
        virtual ~ISurface();

        /**
         * @brief Gets the number of vertices.
         * @return Number of vertices.
         */
        Eigen::Index getNumVertices() const;

        /**
         * @brief Gets the number of dimensions.
         * @return Number of dimensions.
         */
        const Eigen::Index getDimensions() const;

        /**
         * @brief Gets the number of attribute arrays. Returns zero if no attributes have been created yet.
         * @return Number of attributes.
         */
        size_t getNumAttributes() const;

        /**
         * @brief Gets a pointer to the vertex data.
         * @return Array containing vertices.
         */
        std::shared_ptr<IArray> getVertices();

        /**
         * @brief Gets a pointer to the vertex data.
         * @return Array containing vertices.
         */
        std::shared_ptr<const IArray> getVertices() const;

        /**
         * @brief Gets a pointer to the attributes.
         * @return Attribute collection.
         */
        std::shared_ptr<Attributes> getAttributes();

        /**
         * @brief Gets a pointer to the attributes.
         * @return Attribute collection.
         */
        std::shared_ptr<const Attributes> getAttributes() const;

        /**
         * @brief Sets the attributes.
         * @param attributes Attribute collection.
         */
        void setAttributes(std::shared_ptr<Attributes> attributes);

        /**
         * @brief Removes all vertices and attributes.
         */
        virtual void clear();

        /**
         * @brief Recomputes the bounding box from the vertex buffer.
         */
        virtual void recomputeBoundingBox() = 0;

        /**
         * @brief Tests if two surface geometries are equal.
         * @param other Surface to compare with.
         * @return True, if surfaces are equal.
         */
        virtual bool isEqual(const ISurface* other) const;

        /**
         * @brief Gets the number of vertices.
         * @return Number of indices.
         */
        Eigen::Index getNumIndices() const;

        /**
         * @brief Gets a pointer to the index data.
         * @return Index buffer array.
         */
        std::shared_ptr<IArray1> getIndices();

        /**
         * @brief Gets a pointer to the index data.
         * @return Index buffer array.
         */
        std::shared_ptr<const IArray1> getIndices() const;

        /**
         * @brief Gets a pointer to the indices data, which is dynamically casted to a specific indices type. if the requested type does not match the internal type, the function returns null.
         * @tparam TArray Type of index buffer array.
         * @return Index buffer array.
         */
        template <typename TArray>
        std::shared_ptr<TArray> getIndices() { return std::dynamic_pointer_cast<TArray>(mIndices); }

        /**
         * @brief Gets a pointer to the indices data, which is dynamically casted to a specific indices type. if the requested type does not match the internal type, the function returns null.
         * @tparam TArray Type of index buffer array.
         * @return Index buffer array.
         */
        template <typename TArray>
        std::shared_ptr<const TArray> getIndices() const { return std::dynamic_pointer_cast<const TArray>(mIndices); }

        /**
         * @brief Clone the attributes from another point geometry. This deletes all current attributes.
         * @param other Surface to clone indices from.
         */
        void cloneIndices(const ISurface& other);

        /**
         * @brief Removes all indices.
         */
        void clearIndices();

        /**
         * @brief Allocates an array of indices of a specific type. if an array already exists, it is deleted first.
         * @tparam ArrayIndexType Type of index buffer to create.
         * @param length Number of indices.
         * @param primitiveTopology Primitive topology to set.
         * @return Index buffer array.
         */
        template <class ArrayIndexType>
        std::shared_ptr<ArrayIndexType> createIndices(const size_t& length, const ETopology primitiveTopology)
        {
            mPrimitiveTopology = primitiveTopology;
            auto result        = std::make_shared<ArrayIndexType>();
            result->setSize(length);
            mIndices = result;
            return result;
        }

        /**
         * @brief Gets the primitive topology.
         * @return Primitive topology.
         */
        const ETopology& getPrimitiveTopology() const;

        /**
         * @brief Sets the primitive topology.
         * @param topology Primitive topology.
         */
        void setPrimitiveTopology(const ETopology& topology);

    protected:
        /**
         * @brief Gets a pointer to the vertex data.
         * @return Vertex array.
         */
        virtual std::shared_ptr<IArray> getVerticesImpl() = 0;

        /**
         * @brief Gets a pointer to the vertex data.
         * @return Vertex array.
         */
        virtual std::shared_ptr<const IArray> getVerticesImpl() const = 0;

    private:
        /**
         * @brief Attribute that are stored on this geometry.
         */
        std::shared_ptr<Attributes> mAttributes;

        /**
         * @brief Index buffer.
         */
        std::shared_ptr<IArray1> mIndices;

        /**
         * @brief Primitive topology.
         */
        ETopology mPrimitiveTopology;
    };
}
