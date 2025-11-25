#pragma once

#include <vislab/core/data.hpp>
#include <vislab/core/types.hpp>

#include <memory>

namespace vislab
{
    class IArray;
    class Attributes;

    /**
     * @brief Basic interface for point geometry.
     */
    class IPoints : public Data
    {
        VISLAB_INTERFACE(IPoints, Data)

    public:
        /**
         * @brief Constructor.
         */
        IPoints();

        /**
         * @brief Copy-constructor.
         * @param other Points to copy from.
         */
        IPoints(const IPoints& other);

        /**
         * @brief Destructor.
         */
        virtual ~IPoints();

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
         * @brief Appends a point set in the end. If the points have attributes, their types and order must be identical.
         * @param points Point data to append.
         */
        void append(const IPoints* points);

        /**
         * @brief Preprends a point set at the front.
         * @param points Point data to prepend.
         */
        void prepend(const IPoints* points);

        /**
         * @brief Removes the first n elements of this vector.
         * @param n Number of points to remove.
         */
        void removeFirst(size_t n = 1);

        /**
         * @brief Removes the last n elements of this vector.
         * @param n Number of points to remove.
         */
        void removeLast(size_t n = 1);

        /**
         * @brief Reverses the order of the points.
         */
        void reverse();

        /**
         * @brief Recomputes the bounding box from the vertex buffer.
         */
        virtual void recomputeBoundingBox() = 0;

        /**
         * @brief Tests if two point geometries are equal.
         * @param other Points to compare with.
         * @return True, if points are equal.
         */
        virtual bool isEqual(const IPoints* other) const;

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
    };
}
