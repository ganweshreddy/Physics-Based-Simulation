#pragma once

#include "array_fwd.hpp"
#include "base_array.hpp"
#include "traits.hpp"

namespace vislab
{
    /**
     * @brief Type-safe data array with one or multiple components.
     * @tparam TScalar Scalar type for each component.
     * @tparam TDimensions Number of components per element. Use -1 for dynamic number of components.
     * @tparam Dense Dense or sparse array.
     */
    template <typename TScalar, int TDimensions, bool Dense>
    class Array : public BaseArray<TDimensions>
    {
        VISLAB_OBJECT(Array, BaseArray<TDimensions>)

    public:
        /**
         * @brief Underlying Eigen Matrix type.
         */
        using Matrix = std::conditional_t<Dense, Eigen::Matrix<TScalar, TDimensions, -1>, Eigen::SparseMatrix<TScalar>>;

        /**
         * @brief Underlying scalar type.
         */
        using Scalar = TScalar;

        /**
         * @brief Underlying element vector stored at each array index.
         */
        using Element = Eigen::Vector<Scalar, TDimensions>;

        /**
         * @brief Const block reference for returning a const reference to an element.
         */
        using ConstElementRef = typename Matrix::ConstColXpr;

        /**
         * @brief Number of components in each of the element entries.
         */
        static constexpr int Dimensions = TDimensions;

        /**
         * @brief Constructor.
         */
        Array()
        {
            if constexpr (Dimensions != -1)
                mData.conservativeResize(Dimensions, mData.cols());
        }

        /**
         * @brief Gets the underlying Eigen::Matrix. The elements are stored in the columns. The components are stored in the rows. Data is stored in column-major order.
         * @return Eigen::Matrix that stores the data in the array.
         */
        const Matrix& getData() const
        {
            return mData;
        }

        /**
         * @brief Gets the underlying Eigen::Matrix. The elements are stored in the columns. The components are stored in the rows. Data is stored in column-major order.
         * @return Eigen::Matrix that stores the data in the array.
         */
        Matrix& getData()
        {
            return mData;
        }

        /**
         * @brief Gets the number of elements.
         * @return Number of elements in the array.
         */
        const Eigen::Index getSize() const override
        {
            return mData.cols();
        }

        /**
         * @brief Gets the number of components of a single element.
         * @return Number of components per element.
         */
        const Eigen::Index getNumberOfComponents() const override
        {
            return mData.rows();
        }

        /**
         * @brief Get the size in bytes of a single element.
         * @return Size in bytes of each element.
         */
        const Eigen::Index getElementSizeInBytes() const override
        {
            return sizeof(Scalar) * getNumberOfComponents();
        }

        /**
         * @brief Get the size in bytes of the entire array (only counting the data, not meta information like the array name).
         * @return Total size in bytes.
         */
        const Eigen::Index getSizeInBytes() const override
        {
            return getSize() * getElementSizeInBytes();
        }

        /**
         * @brief Sets the number of elements.
         * @param size New number of elements.
         */
        void setSize(const Eigen::Index& size) override
        {
            mData.conservativeResize(mData.rows(), size);
        }

        /**
         * @brief Sets the number of components. This function is only available for arrays with variable number of components.
         * @tparam D SFINAE parameter.
         * @param numComponents New number of components.
         * @return void
         */
        template <int D = Dimensions>
        typename std::enable_if_t<D == -1>
        setNumberOfComponents(const Eigen::Index& numComponents)
        {
            mData.conservativeResize(numComponents, mData.cols());
        }

        /**
         * @brief Gets the element at a specific index in its internal format.
         * @param index Array index.
         * @return Element at a certain array index.
         */
        ConstElementRef getValue(const Eigen::Index& index) const
        {
            return mData.col(index);
        }

        /**
         * @brief Sets the element at a specific index with all its components in its internal format.
         * @param index Array index.
         * @param value New element to set at the array index.
         */
        void setValue(const Eigen::Index& index, const Element& value)
        {
            assert(mData.rows() == value.rows()); // make sure that the element "value" has the right dimension
            if constexpr (detail::is_sparse_eigen_matrix<Matrix>())
                mData.col(index) = value.sparseView();
            else
                mData.col(index) = value;
        }

        /**
         * @brief Sets the element at a specific index if it is a scalar value.
         * @tparam D SFINAE parameter.
         * @param index Array index.
         * @param value New scalar element to set at the array index.
         * @return void
         */
        template <int D = Dimensions>
        typename std::enable_if_t<(D == 1)> setValue(const Eigen::Index& index, const Scalar& value)
        {
            mData.coeffRef(0, index) = value;
        }

        /**
         * @brief Gets a value in double precision by casting from the internal format. Note that this access involves a type conversion, which entails a copy.
         * @param index Array index.
         * @return Element in double format.
         */
        Eigen::Vector<double, Dimensions> getValueDouble(const Eigen::Index& index) const
        {
            return getValueDoubleImplStatic(index);
        }

        /**
         * @brief Sets a value in double precision by casting to the internal format.
         * @param index Array index.
         * @param value Element to set in double format.
         */
        void setValueDouble(const Eigen::Index& index, const Eigen::Vector<double, Dimensions>& value)
        {
            setValueDoubleImplStatic(index, value);
        }

        /**
         * @brief Sets the element at a specific index in double precision by casting to the internal format if it is a scalar value.
         * @tparam D SFINAE parameter.
         * @param index Array index.
         * @param value New scalar element to set at the array index.
         * @return void
         */
        template <int D = Dimensions>
        typename std::enable_if_t<(D == 1)> setValueDouble(const Eigen::Index& index, const double& value)
        {
            mData.coeffRef(0, index) = static_cast<Scalar>(value);
        }

        /**
         * @brief Gets a value in int64_t precision by casting from the internal format. Note that this access involves a type conversion, which entails a copy.
         * @param index Array index.
         * @return Element in int64_t format.
         */
        Eigen::Vector<int64_t, Dimensions> getValueInt64(const Eigen::Index& index) const
        {
            return getValueInt64ImplStatic(index);
        }

        /**
         * @brief Sets a value in int64_t precision by casting to the internal format.
         * @param index Array index.
         * @param value Element to set in int64_t format.
         */
        void setValueInt64(const Eigen::Index& index, const Eigen::Vector<int64_t, Dimensions>& value)
        {
            setValueInt64ImplStatic(index, value);
        }

        /**
         * @brief Sets the element at a specific index in double precision by casting to the internal format if it is a scalar value.
         * @tparam D SFINAE parameter.
         * @param index Array index.
         * @param value New scalar element to set at the array index.
         * @return void
         */
        template <int D = Dimensions>
        typename std::enable_if_t<(D == 1)> setValueInt64(const Eigen::Index& index, const int64_t& value)
        {
            mData.coeffRef(0, index) = static_cast<Scalar>(value);
        }

        /**
         * @brief Constructs the array with an initializer list.
         * @tparam D SFINAE parameter.
         * @param list List of values.
         * @return void
         */
        template <int D = Dimensions>
        typename std::enable_if_t<(D == 1)> setValues(std::initializer_list<TScalar> list)
        {
            setSize(list.size());
            size_t i = 0;
            for (auto elem : list)
                setValue(i++, elem);
        }

        /**
         * @brief Constructs the array with an initializer list.
         * @param list List of elements.
         */
        void setValues(std::initializer_list<Element> list)
        {
            setSize(list.size());
            size_t i = 0;
            for (auto elem : list)
                setValue(i++, elem);
        }

        /**
         * @brief Sets all values to zero.
         */
        void setZero() override
        {
            mData.setZero();
        }

        /**
         * @brief Deletes all elements.
         */
        void clear() override
        {
            setSize(0);
        }

        /**
         * @brief Appends a value at the end.
         * @param value Value to append.
         */
        void append(const Element& value)
        {
            setSize(getSize() + 1);
            setValue(getSize() - 1, value);
        }

        /**
         * @brief Appends an entire array at the end.
         * @param other Array to append.
         */
        void append(const Array& other)
        {
            size_t oldsize = getSize();
            setSize(oldsize + other.getSize());
            for (int i = 0; i < other.getSize(); ++i)
                setValue(oldsize + i, other.getValue(i));
        }

        /**
         * @brief Appends an entire array at the end. Must be the same type, otherwise this function does nothing.
         * @param other Array to append.
         */
        void append(const IArray* other) override
        {
            if (dynamic_cast<const Array*>(other))
                append(*dynamic_cast<const Array*>(other));
        }

        /**
         * @brief Prepends an element at the front (linear time complexity, since the entire array is copied).
         * @param value Element to prepend.
         */
        void prepend(const Element& value)
        {
            size_t oldsize = getSize();
            setSize(getSize() + 1);
            for (int i = (int)oldsize - 1; i >= 0; --i)
                setValue(i + 1, getValue(i));
            setValue(0, value);
        }

        /**
         * @brief Prepends an entire array at the front (linear time complexity, since the entire array is copied).
         * @param other Array to prepend.
         */
        void prepend(const Array& other)
        {
            // resize the buffer to accomodate the current and the other ones data
            size_t oldsize   = getSize();
            size_t othersize = other.getSize();
            setSize(oldsize + othersize);
            // move all entries to the end
            for (int i = (int)oldsize - 1; i >= 0; --i)
                setValue(othersize + i, getValue(i));
            // insert the new values at the front
            for (size_t i = 0; i < othersize; ++i)
                setValue(i, other.getValue(i));
        }

        /**
         * @brief Prepends an entire array at the front (linear time complexity, since the entire array is copied). Must be the same type, otherwise this function does nothing.
         * @param other Array to prepend.
         */
        void prepend(const IArray* other) override
        {
            if (dynamic_cast<const Array*>(other))
                prepend(*dynamic_cast<const Array*>(other));
        }

        /**
         * @brief Reverses the order of the elements.
         */
        void reverse() override
        {
            const Eigen::Index n = getSize();
            for (Eigen::Index i = 0; i < n / 2; ++i)
            {
                Element temp = getValue(i);
                setValue(i, getValue(n - 1 - i));
                setValue(n - 1 - i, temp);
            }
        }

        /**
         * @brief Removes the last n elements of this vector.
         * @param n Number of elements to remove at the end.
         */
        void removeLast(size_t n = 1) override
        {
            setSize(getSize() - n);
        }

        /**
         * @brief Removes the first n elements of this vector.
         * @param n Number of elements to remove at the front.
         */
        void removeFirst(size_t n = 1) override
        {
            // Copy each element n places ahead
            for (size_t i = n; i < getSize(); ++i)
                setValue(i - n, getValue(i));
            removeLast(n);
        }

        /**
         * @brief Get the first element.
         * @return First element.
         */
        ConstElementRef first() const
        {
            return getValue(0);
        }

        /**
         * @brief Get the last element.
         * @return Last element.
         */
        ConstElementRef last() const
        {
            return getValue(getSize() - 1);
        }

        /**
         * @brief Gets component-wise the minimum value (computed in linear time).
         * @return Smallest value per component.
         */
        Element getMin() const
        {
            Element val = getValue(0);
            for (int i = 0; i < Dimensions; ++i)
            {
                if constexpr (detail::is_sparse_eigen_matrix<Matrix>())
                {
                    TScalar minValue = std::numeric_limits<TScalar>::max();
                    for (typename vislab::Array<TScalar, Dimensions, Dense>::Matrix::InnerIterator it(mData, i); it; ++it)
                        minValue = std::min(minValue, it.value());
                    val(i) = minValue;
                }
                else
                    val(i) = mData.row(i).minCoeff(&val[i]);
            }
            return val;
        }

        /**
         * @brief Gets component-wise the maximum value (computed in linear time).
         * @return Largest value per component.
         */
        Element getMax() const
        {
            Element val = getValue(0);
            for (int i = 0; i < Dimensions; ++i)
            {
                if constexpr (detail::is_sparse_eigen_matrix<Matrix>())
                {
                    TScalar maxValue = -std::numeric_limits<TScalar>::max();
                    for (typename vislab::Array<TScalar, Dimensions, Dense>::Matrix::InnerIterator it(mData, i); it; ++it)
                        maxValue = std::max(maxValue, it.value());
                    val(i) = maxValue;
                }
                else
                    val(i) = mData.row(i).maxCoeff(&val[i]);
            }
            return val;
        }

        /**
         * @brief Sort componentwise in ascending order.
         */
        void sortAscending() override
        {
            std::vector<Scalar> entries(getSize());
            for (int c = 0; c < getNumberOfComponents(); ++c)
            {
                for (Eigen::Index i = 0; i < getSize(); ++i)
                    entries[i] = mData.coeffRef(c, i);
                std::sort(entries.begin(), entries.end());
                for (Eigen::Index i = 0; i < getSize(); ++i)
                    mData.coeffRef(c, i) = entries[i];
            }
        }

        /**
         * @brief Sort componentwise in descending order.
         */
        void sortDescending() override
        {
            std::vector<Scalar> entries(getSize());
            for (int c = 0; c < getNumberOfComponents(); ++c)
            {
                for (Eigen::Index i = 0; i < getSize(); ++i)
                    entries[i] = mData.coeffRef(c, i);
                std::sort(entries.rbegin(), entries.rend());
                for (Eigen::Index i = 0; i < getSize(); ++i)
                    mData.coeffRef(c, i) = entries[i];
            }
        }

        /**
         * @brief Tests if two arrays are equal.
         * @param other Array to compare with.
         * @return True if equal.
         */
        bool isEqual(const IArray* other) const override
        {
            const Array* otherTyped = dynamic_cast<const Array*>(other);
            if (!otherTyped)
                return false;
            if (this->getSize() != otherTyped->getSize())
                return false;
            if (this->getName() != other->getName())
                return false;
            return mData.isApprox(otherTyped->mData);
        }

        /**
         * @brief Conversion into an Eigen::MatrixXd. You can also use GetData() to get the native type and handle the conversion yourself.
         * @param output Eigen matrix with variable size.
         */
        void convertTo(Eigen::MatrixXd& output) const override
        {
            output = mData.template cast<double>();
        }

        /**
         * @brief Conversion from an Eigen::MatrixXd. You can also use GetData() to get the native type and handle the conversion yourself.
         * @param input Eigen matrix with variable size.
         */
        void convertFrom(const Eigen::MatrixXd& input) override
        {
            if constexpr (detail::is_sparse_eigen_matrix<Matrix>())
                mData = input.template cast<Scalar>().sparseView();
            else
                mData = input.template cast<Scalar>();
        }

    protected:
        /**
         * @brief Gets an element in double precision by casting from the internal format.
         * @param index Array index.
         * @return Element in double precision.
         */
        Eigen::Vector<double, Dimensions> getValueDoubleImplStatic(const Eigen::Index& index) const
        {
            if constexpr (std::is_same_v<Scalar, double>)
                return mData.col(index);
            else
                return mData.col(index).template cast<double>(); // numeric_cast<Eigen::Vector<double, Dimensions>, const Eigen::Vector<Scalar, Dimensions>>();
        }

        /**
         * @brief Sets an element in double precision by casting to the internal format.
         * @param index Array index.
         * @param value Element in double precision.
         */
        void setValueDoubleImplStatic(const Eigen::Index& index, const Eigen::Vector<double, Dimensions>& value)
        {
            assert(mData.rows() == value.rows()); // make sure that the element "value" has the right dimension
            if constexpr (detail::is_sparse_eigen_matrix<Matrix>())
                mData.col(index) = value.template cast<Scalar>().sparseView();
            else
                mData.col(index) = value.template cast<Scalar>();
        }

        /**
         * @brief Gets an element in int64_t precision by casting from the internal format.
         * @param index Array index.
         * @return Element in int64_t precision.
         */
        Eigen::Vector<int64_t, Dimensions> getValueInt64ImplStatic(const Eigen::Index& index) const
        {
            return mData.col(index).template cast<int64_t>();
        }

        /**
         * @brief Sets an element in int64_t precision by casting to the internal format.
         * @param index Array index.
         * @param value Element in int64_t precision.
         */
        void setValueInt64ImplStatic(const Eigen::Index& index, const Eigen::Vector<int64_t, Dimensions>& value)
        {
            assert(mData.rows() == value.rows()); // make sure that the element "value" has the right dimension
            if constexpr (detail::is_sparse_eigen_matrix<Matrix>())
                mData.col(index) = value.template cast<Scalar>().sparseView();
            else
                mData.col(index) = value.template cast<Scalar>();
        }

        /**
         * @brief Gets an element in double precision by casting from the internal format.
         * @param index Array index.
         * @return Element in double precision.
         */
        Eigen::VectorXd getValueDoubleImpl(const Eigen::Index& index) const override
        {
            return getValueDoubleImplStatic(index);
        }

        /**
         * @brief Sets an element in double precision by casting to the internal format.
         * @param index Array index.
         * @param value Element in double precision.
         */
        void setValueDoubleImpl(const Eigen::Index& index, const Eigen::VectorXd& value) override
        {
            setValueDoubleImplStatic(index, value);
        }

        /**
         * @brief Gets a value in int64_t precision by casting from the internal format.
         * @param index Array index.
         * @return Element in int64_t precision.
         */
        Eigen::VectorXl getValueInt64Impl(const Eigen::Index& index) const override
        {
            return getValueInt64ImplStatic(index);
        }

        /**
         * @brief Sets a value in int64_t precision by casting to the internal format.
         * @param index Array index.
         * @param value Element in int64_t precision.
         */
        void setValueInt64Impl(const Eigen::Index& index, const Eigen::VectorXl& value) override
        {
            setValueInt64ImplStatic(index, value);
        }

        /**
         * @brief Underlying Eigen matrix containing the data.
         */
        Matrix mData;
    };
}
