#pragma once

#include "data.hpp"
#include "types.hpp"

namespace vislab
{
    /**
     * @brief Abstract base class for data arrays.
     */
    class IArray : public Data
    {
        VISLAB_INTERFACE(IArray, Data)

    public:
        /**
         * @brief Constructor for an abstract data array.
         */
        IArray();

        /**
         * @brief Gets the name of this data array.
         * @return Name of the data array.
         */
        const std::string& getName() const;

        /**
         * @brief Sets the name of this data array.
         * @param name Name of the data array.
         */
        void setName(const std::string& name);

        /**
         * @brief Gets the number of elements.
         * @return Number of elements in the array.
         */
        virtual const Eigen::Index getSize() const = 0;

        /**
         * @brief Gets the number of components of a single element.
         * @return Number of components per element.
         */
        virtual const Eigen::Index getNumberOfComponents() const = 0;

        /**
         * @brief Get the size in bytes of a single element.
         * @return Size in bytes of each element.
         */
        virtual const Eigen::Index getElementSizeInBytes() const = 0;

        /**
         * @brief Get the size in bytes of the entire array (only counting the data, not meta information like the array name).
         * @return Total size in bytes.
         */
        virtual const Eigen::Index getSizeInBytes() const = 0;

        /**
         * @brief Sets the number of elements.
         * @param size New number of elements.
         */
        virtual void setSize(const Eigen::Index& newSize) = 0;

        /**
         * @brief Deletes all elements.
         */
        virtual void clear() = 0;

        /**
         * @brief Appends an entire array at the end. Must be the same type, otherwise this function does nothing.
         * @param other Array to append.
         */
        virtual void append(const IArray* other) = 0;

        /**
         * @brief Prepends an entire array at the front (linear time complexity, since the entire array is copied). Must be the same type, otherwise this function does nothing.
         * @param other Array to prepend.
         */
        virtual void prepend(const IArray* other) = 0;

        /**
         * @brief Reverses the order of the elements.
         */
        virtual void reverse() = 0;

        /**
         * @brief Removes the last n elements of this vector.
         * @param n Number of elements to remove at the end.
         */
        virtual void removeLast(size_t n = 1) = 0;

        /**
         * @brief Removes the first n elements of this vector.
         * @param n Number of elements to remove at the front.
         */
        virtual void removeFirst(size_t n = 1) = 0;

        /**
         * @brief Sort componentwise in ascending order.
         */
        virtual void sortAscending() = 0;

        /**
         * @brief Sort componentwise in descending order.
         */
        virtual void sortDescending() = 0;

        /**
         * @brief Tests if two arrays are equal.
         * @param other Array to compare with.
         * @return True if equal.
         */
        virtual bool isEqual(const IArray* other) const = 0;

        /**
         * @brief Sets all values to zero.
         */
        virtual void setZero() = 0;

        /**
         * @brief Gets a value in double precision by casting from the internal format. Note that this access involves a type conversion, which entails a copy.
         * @param index Array index.
         * @return Element in double format.
         */
        Eigen::VectorXd getValueDouble(const Eigen::Index& index) const;

        /**
         * @brief Sets a value in double precision by casting to the internal format.
         * @param index Array index.
         * @param value Element to set in double format.
         */
        void setValueDouble(const Eigen::Index& index, const Eigen::VectorXd& value);

        /**
         * @brief Gets a value in int64_t precision by casting from the internal format. Note that this access involves a type conversion, which entails a copy.
         * @param index Array index.
         * @return Element in int64_t format.
         */
        Eigen::VectorXl getValueInt64(const Eigen::Index& index) const;

        /**
         * @brief Sets a value in int64_t precision by casting to the internal format.
         * @param index Array index.
         * @param value Element to set in int64_t format.
         */
        void setValueInt64(const Eigen::Index& index, const Eigen::VectorXl& value);

        /**
         * @brief Conversion into an Eigen::MatrixXd. You can also use GetData() to get the native type and handle the conversion yourself.
         * @param output Eigen matrix with variable size.
         */
        virtual void convertTo(Eigen::MatrixXd& output) const = 0;

        /**
         * @brief Conversion from an Eigen::MatrixXd. You can also use GetData() to get the native type and handle the conversion yourself.
         * @param input Eigen matrix with variable size.
         */
        virtual void convertFrom(const Eigen::MatrixXd& input) = 0;

    protected:
        /**
         * @brief Gets an element in double precision by casting from the internal format.
         * @param index Array index.
         * @return Element in double precision.
         */
        virtual Eigen::VectorXd getValueDoubleImpl(const Eigen::Index& index) const = 0;

        /**
         * @brief Sets an element in double precision by casting to the internal format.
         * @param index Array index.
         * @param value Element in double precision.
         */
        virtual void setValueDoubleImpl(const Eigen::Index& index, const Eigen::VectorXd& value) = 0;

        /**
         * @brief Gets a value in int64_t precision by casting from the internal format.
         * @param index Array index.
         * @return Element in int64_t precision.
         */
        virtual Eigen::VectorXl getValueInt64Impl(const Eigen::Index& index) const = 0;

        /**
         * @brief Sets a value in int64_t precision by casting to the internal format.
         * @param index Array index.
         * @param value Element in int64_t precision.
         */
        virtual void setValueInt64Impl(const Eigen::Index& index, const Eigen::VectorXl& value) = 0;

    private:
        /**
         * @brief Name of the array.
         */
        std::string mName;
    };
}
