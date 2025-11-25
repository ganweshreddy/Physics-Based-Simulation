#pragma once

#include "lines_fwd.hpp"

#include "base_lines.hpp"
#include "line.hpp"

#include <vislab/core/array.hpp>
#include <vislab/core/traits.hpp>

namespace vislab
{
    /**
     * @brief Set of line geometry with specific internal type.
     * @tparam TArrayType Internal array type used to store vertices.
     */
    template <typename TArrayType>
    class Lines : public BaseLines<TArrayType::Dimensions>
    {
        VISLAB_OBJECT(Lines, BaseLines<TArrayType::Dimensions>)

    public:
        /**
         * @brief Type of the bounding box.
         */
        using BoundingBox = Eigen::AlignedBox<double, TArrayType::Dimensions>;

        /**
         * @brief Number of dimensions for vertex positions.
         */
        static const int Dimensions = TArrayType::Dimensions;

        /**
         * @brief Internal array type used to store the vertices.
         */
        using ArrayType = TArrayType;

        /**
         * @brief Internal line type.
         */
        using LineType = Line<TArrayType>;

        /**
         * @brief Constructor.
         */
        Lines()
        {
        }

        /**
         * @brief Copy-constructor.
         * @param other Geometry to copy from.
         */
        Lines(const Lines& other)
            : BaseLines<TArrayType::Dimensions>(other)
        {
            for (auto line : other.lines)
                lines.push_back(std::shared_ptr<LineType>(line->clone()));
        }

        /**
         * @brief Destructor.
         */
        virtual ~Lines() { clear(); }

        /**
         * @brief Creates a new line and returns a pointer to it. This function is thread-safe.
         * @return Created line.
         */
        std::shared_ptr<LineType> createLine()
        {
            std::shared_ptr<LineType> line;
#ifdef NDEBUG
#pragma omp critical
#endif
            {
                lines.emplace_back(std::make_shared<LineType>());
                line = lines.back();
            }
            return line;
        }

        /**
         * @brief Creates a new line, which is copied from an existing line and returns a pointer to it. This function is thread-safe.
         * @param copy Line to make a copy from.
         * @return Created line.
         */
        std::shared_ptr<LineType> createLine(const LineType& copy)
        {
            std::shared_ptr<LineType> line;
#ifdef NDEBUG
#pragma omp critical
#endif
            {
                lines.emplace_back(std::make_shared<LineType>(copy));
                line = lines.back();
            }
            return line;
        }

        /**
         * @brief Gets the number of lines in this container.
         * @return Number of lines.
         */
        size_t getNumberOfLines() const override { return lines.size(); }

        /**
         * @brief Deletes all lines
         */
        void clear() override
        {
            lines.clear();
        }

        /**
         * @brief Gets a line.
         * @param index Index of line to get.
         * @return Line to get.
         */
        std::shared_ptr<LineType> getLine(size_t index)
        {
            return std::dynamic_pointer_cast<LineType>(getLineImpl(index));
        }

        /**
         * @brief Gets a line.
         * @param index Index of line to get.
         * @return Line to get.
         */
        std::shared_ptr<const LineType> getLine(size_t index) const
        {
            return std::dynamic_pointer_cast<const LineType>(getLineImpl(index));
        }

        /**
         * @brief Adds a line that matches the type of the derived class.
         * @param line Line to add.
         */
        void addLine(std::shared_ptr<LineType> line)
        {
#ifdef NDEBUG
#pragma omp critical
#endif
            {
                lines.push_back(line);
            }
        }

        /**
         * @brief Collection of all the lines in the line set.
         */
        std::vector<std::shared_ptr<LineType>> lines;

        /**
         * @brief Tests if two line geometries are equal.
         * @param other Line set to compare with.
         * @return True, if the line sets are equal.
         */
        bool isEqual(const ILines* other) const override
        {
            if (!BaseLines<TArrayType::Dimensions>::isEqual(other))
                return false;
            auto otherTyped = dynamic_cast<const Lines*>(other);
            if (otherTyped == nullptr)
                return false;
            if (lines.size() != otherTyped->lines.size())
                return false;
            for (size_t i = 0; i < lines.size(); ++i)
                if (!lines[i]->isEqual(otherTyped->lines[i].get()))
                    return false;
            return true;
        }

    private:
        /**
         * @brief Gets a line.
         * @param index Index of the line to get.
         * @return Line to get.
         */
        std::shared_ptr<ILine> getLineImpl(size_t index) override { return lines[index]; }

        /**
         * @brief Gets a line.
         * @param index Index of the line to get.
         * @return Line to get.
         */
        std::shared_ptr<const ILine> getLineImpl(size_t index) const override { return lines[index]; }

        /**
         * @brief Adds a line that matches the LineGeometry::Line type of the derived class.
         * @param line Line to add.
         */
        void addLineImpl(std::shared_ptr<ILine> line) override
        {
            auto l = std::dynamic_pointer_cast<LineType>(line);
            if (l != nullptr)
                lines.push_back(l);
        }
    };
}
