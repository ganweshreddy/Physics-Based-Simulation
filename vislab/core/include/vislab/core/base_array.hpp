#pragma once

#include "base_array_fwd.hpp"

#include "iarray.hpp"

namespace vislab
{
    /**
     * @brief Abstract base class for data arrays with a fixed number of components.
     * @tparam TDimensions Number of components per element. Use -1 for dynamic number of components.
     */
    template <int64_t TDimensions>
    class BaseArray : public IArray
    {
        VISLAB_INTERFACE(BaseArray, IArray)

    public:
        /**
         * @brief Number of components in each of the element entries.
         */
        static const int64_t Dimensions = TDimensions;

        /**
         * @brief Gets an element in double precision by casting from the internal format. Note that this access involves a type conversion, which entails a copy.
         * @param index Array index.
         * @return Element in double precision.
         */
        Eigen::Vector<double, TDimensions> getValueDouble(const Eigen::Index& index) const
        {
            return getValueDoubleImpl(index);
        }

        /**
         * @brief Sets a value in double precision by casting to the internal format.
         * @param index Array index.
         * @param value Element in double precision.
         */
        void setValueDouble(const Eigen::Index& index, const Eigen::Vector<double, TDimensions>& value)
        {
            setValueDoubleImpl(index, value);
        }

        /**
         * @brief Gets a value in int64_t precision by casting from the internal format.  Note that this access involves a type conversion, which entails a copy.
         * @param index Array index.
         * @return Element in int64_t precision.
         */
        Eigen::Vector<int64_t, TDimensions> getValueInt64(const Eigen::Index& index) const
        {
            return getValueInt64Impl(index);
        }

        /**
         * @brief Sets a value in int64_t precision by casting to the internal format.
         * @param index Array index.
         * @param value Element in int64_t precision.
         */
        void setValueInt64(const Eigen::Index& index, const Eigen::Vector<int64_t, TDimensions>& value)
        {
            setValueInt64Impl(index, value);
        }
    };
}
