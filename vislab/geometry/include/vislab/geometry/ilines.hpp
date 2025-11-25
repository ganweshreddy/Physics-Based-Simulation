#pragma once

#include "iline.hpp"

namespace vislab
{
    /**
     * @brief Basic interface for set of line geometry.
     */
    class ILines : public Data
    {
        VISLAB_INTERFACE(ILines, Data)

    public:
        /**
         * @brief Gets the number of lines in this container.
         * @return Number of lines.
         */
        virtual size_t getNumberOfLines() const = 0;

        /**
         * @brief Deletes all lines (calls destructors)
         */
        virtual void clear() = 0;

        /**
         * @brief Recomputes the bounding box from the vertex buffers.
         */
        virtual void recomputeBoundingBox() = 0;

        /**
         * @brief Gets a line.
         * @param index Index of the line to get.
         * @return Requested line.
         */
        std::shared_ptr<ILine> getLine(size_t index);

        /**
         * @brief Gets a line.
         * @param index Index of the line to get.
         * @return Requested line.
         */
        std::shared_ptr<const ILine> getLine(size_t index) const;

        /**
         * @brief Adds a line that matches the type of the derived class.
         * @param line Line to add.
         */
        void addLine(std::shared_ptr<ILine> line);

        /**
         * @brief Tests if two line geometries are equal.
         * @param other Line set to compare with.
         * @return True, if line sets are equal.
         */
        virtual bool isEqual(const ILines* other) const = 0;

    protected:
        /**
         * @brief Gets a line.
         * @param index Index of the line to get.
         * @return Requested line.
         */
        virtual std::shared_ptr<ILine> getLineImpl(size_t index) = 0;

        /**
         * @brief Gets a line.
         * @param index Index of the line to get.
         * @return Requested line.
         */
        virtual std::shared_ptr<const ILine> getLineImpl(size_t index) const = 0;

        /**
         * @brief Adds a line that matches the type of the derived class.
         * @param line Line to add.
         */
        virtual void addLineImpl(std::shared_ptr<ILine> line) = 0;
    };
}
