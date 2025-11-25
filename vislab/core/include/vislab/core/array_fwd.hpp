#pragma once

#include <stdint.h>

namespace vislab
{
    /**
     * @brief Type-safe data array with one or multiple components. Use Dimensions=-1 for dynamic number of components.
     * @tparam TScalar Scalar type that is stored in the array.
     * @tparam TDimensions Number of components per element.
     * @tparam Dense Flag that specifies whether this is a dense or sparse array.
     */
    template <typename TScalar, int TDimensions, bool Dense = true>
    class Array;

    /**
     * @brief Dense one-dimensional array of float values with one component.
     */
    using Array1f = Array<float, 1, true>;

    /**
     * @brief Dense one-dimensional array of float values with two components.
     */
    using Array2f = Array<float, 2, true>;

    /**
     * @brief Dense one-dimensional array of float values with three components.
     */
    using Array3f = Array<float, 3, true>;

    /**
     * @brief Dense one-dimensional array of float values with four components.
     */
    using Array4f = Array<float, 4, true>;

    /**
     * @brief Dense one-dimensional array of float values with variable number of components.
     */
    using ArrayXf = Array<float, -1, true>;

    /**
     * @brief Dense one-dimensional array of double values with one component.
     */
    using Array1d = Array<double, 1, true>;

    /**
     * @brief Dense one-dimensional array of double values with two components.
     */
    using Array2d = Array<double, 2, true>;

    /**
     * @brief Dense one-dimensional array of double values with three components.
     */
    using Array3d = Array<double, 3, true>;

    /**
     * @brief Dense one-dimensional array of double values with four components.
     */
    using Array4d = Array<double, 4, true>;

    /**
     * @brief Dense one-dimensional array of double values with variable number of components.
     */
    using ArrayXd = Array<double, -1, true>;

    /**
     * @brief Dense one-dimensional array of int16_t values with one component.
     */
    using Array1i_16 = Array<int16_t, 1, true>;

    /**
     * @brief Dense one-dimensional array of int16_t values with two components.
     */
    using Array2i_16 = Array<int16_t, 2, true>;

    /**
     * @brief Dense one-dimensional array of int16_t values with three components.
     */
    using Array3i_16 = Array<int16_t, 3, true>;

    /**
     * @brief Dense one-dimensional array of int16_t values with four components.
     */
    using Array4i_16 = Array<int16_t, 4, true>;

    /**
     * @brief Dense one-dimensional array of int16_t values with variable number of components.
     */
    using ArrayXi_16 = Array<int16_t, -1, true>;

    /**
     * @brief Dense one-dimensional array of int32_t values with one component.
     */
    using Array1i = Array<int32_t, 1, true>;

    /**
     * @brief Dense one-dimensional array of int32_t values with two components.
     */
    using Array2i = Array<int32_t, 2, true>;

    /**
     * @brief Dense one-dimensional array of int32_t values with three components.
     */
    using Array3i = Array<int32_t, 3, true>;

    /**
     * @brief Dense one-dimensional array of int32_t values with four components.
     */
    using Array4i = Array<int32_t, 4, true>;

    /**
     * @brief Dense one-dimensional array of int32_t values with variable number of components.
     */
    using ArrayXi = Array<int32_t, -1, true>;

    /**
     * @brief Dense one-dimensional array of int64_t values with one component.
     */
    using Array1i_64 = Array<int64_t, 1, true>;

    /**
     * @brief Dense one-dimensional array of int64_t values with two components.
     */
    using Array2i_64 = Array<int64_t, 2, true>;

    /**
     * @brief Dense one-dimensional array of int64_t values with three components.
     */
    using Array3i_64 = Array<int64_t, 3, true>;

    /**
     * @brief Dense one-dimensional array of int64_t values with four components.
     */
    using Array4i_64 = Array<int64_t, 4, true>;

    /**
     * @brief Dense one-dimensional array of int64_t values with variable number of components.
     */
    using ArrayXi_64 = Array<int64_t, -1, true>;

    /**
     * @brief Dense one-dimensional array of uint16_t values with one component.
     */
    using Array1ui_16 = Array<uint16_t, 1, true>;

    /**
     * @brief Dense one-dimensional array of uint16_t values with two components.
     */
    using Array2ui_16 = Array<uint16_t, 2, true>;

    /**
     * @brief Dense one-dimensional array of uint16_t values with three components.
     */
    using Array3ui_16 = Array<uint16_t, 3, true>;

    /**
     * @brief Dense one-dimensional array of uint16_t values with four components.
     */
    using Array4ui_16 = Array<uint16_t, 4, true>;

    /**
     * @brief Dense one-dimensional array of uint16_t values with variable number of components.
     */
    using ArrayXui_16 = Array<uint16_t, -1, true>;

    /**
     * @brief Dense one-dimensional array of uint32_t values with one component.
     */
    using Array1ui = Array<uint32_t, 1, true>;

    /**
     * @brief Dense one-dimensional array of uint32_t values with two components.
     */
    using Array2ui = Array<uint32_t, 2, true>;

    /**
     * @brief Dense one-dimensional array of uint32_t values with three components.
     */
    using Array3ui = Array<uint32_t, 3, true>;

    /**
     * @brief Dense one-dimensional array of uint32_t values with four components.
     */
    using Array4ui = Array<uint32_t, 4, true>;

    /**
     * @brief Dense one-dimensional array of uint32_t values with variable number of components.
     */
    using ArrayXui = Array<uint32_t, -1, true>;

    /**
     * @brief Dense one-dimensional array of uint64_t values with one component.
     */
    using Array1ui_64 = Array<uint64_t, 1, true>;

    /**
     * @brief Dense one-dimensional array of uint64_t values with two components.
     */
    using Array2ui_64 = Array<uint64_t, 2, true>;

    /**
     * @brief Dense one-dimensional array of uint64_t values with three components.
     */
    using Array3ui_64 = Array<uint64_t, 3, true>;

    /**
     * @brief Dense one-dimensional array of uint64_t values with four components.
     */
    using Array4ui_64 = Array<uint64_t, 4, true>;

    /**
     * @brief Dense one-dimensional array of uint64_t values with variable number of components.
     */
    using ArrayXui_64 = Array<uint64_t, -1, true>;

    /**
     * @brief Sparse one-dimensional array of float values with one component.
     */
    using SparseArray1f = Array<float, 1, false>;

    /**
     * @brief Sparse one-dimensional array of float values with two components.
     */
    using SparseArray2f = Array<float, 2, false>;

    /**
     * @brief Sparse one-dimensional array of float values with three components.
     */
    using SparseArray3f = Array<float, 3, false>;

    /**
     * @brief Sparse one-dimensional array of float values with four components.
     */
    using SparseArray4f = Array<float, 4, false>;

    /**
     * @brief Sparse one-dimensional array of float values with variable number of components.
     */
    using SparseArrayXf = Array<float, -1, false>;

    /**
     * @brief Sparse one-dimensional array of double values with one component.
     */
    using SparseArray1d = Array<double, 1, false>;

    /**
     * @brief Sparse one-dimensional array of double values with two components.
     */
    using SparseArray2d = Array<double, 2, false>;

    /**
     * @brief Sparse one-dimensional array of double values with three components.
     */
    using SparseArray3d = Array<double, 3, false>;

    /**
     * @brief Sparse one-dimensional array of double values with four components.
     */
    using SparseArray4d = Array<double, 4, false>;

    /**
     * @brief Sparse one-dimensional array of double values with variable number of components.
     */
    using SparseArrayXd = Array<double, -1, false>;

    /**
     * @brief Sparse one-dimensional array of int16_t values with one component.
     */
    using SparseArray1i_16 = Array<int16_t, 1, false>;

    /**
     * @brief Sparse one-dimensional array of int16_t values with two components.
     */
    using SparseArray2i_16 = Array<int16_t, 2, false>;

    /**
     * @brief Sparse one-dimensional array of int16_t values with three components.
     */
    using SparseArray3i_16 = Array<int16_t, 3, false>;

    /**
     * @brief Sparse one-dimensional array of int16_t values with four components.
     */
    using SparseArray4i_16 = Array<int16_t, 4, false>;

    /**
     * @brief Sparse one-dimensional array of int16_t values with variable number of components.
     */
    using SparseArrayXi_16 = Array<int16_t, -1, false>;

    /**
     * @brief Sparse one-dimensional array of int32_t values with one component.
     */
    using SparseArray1i = Array<int32_t, 1, false>;

    /**
     * @brief Sparse one-dimensional array of int32_t values with two components.
     */
    using SparseArray2i = Array<int32_t, 2, false>;

    /**
     * @brief Sparse one-dimensional array of int32_t values with three components.
     */
    using SparseArray3i = Array<int32_t, 3, false>;

    /**
     * @brief Sparse one-dimensional array of int32_t values with four components.
     */
    using SparseArray4i = Array<int32_t, 4, false>;

    /**
     * @brief Sparse one-dimensional array of int32_t values with variable number of components.
     */
    using SparseArrayXi = Array<int32_t, -1, false>;

    /**
     * @brief Sparse one-dimensional array of int64_t values with one component.
     */
    using SparseArray1i_64 = Array<int64_t, 1, false>;

    /**
     * @brief Sparse one-dimensional array of int64_t values with two components.
     */
    using SparseArray2i_64 = Array<int64_t, 2, false>;

    /**
     * @brief Sparse one-dimensional array of int64_t values with three components.
     */
    using SparseArray3i_64 = Array<int64_t, 3, false>;

    /**
     * @brief Sparse one-dimensional array of int64_t values with four components.
     */
    using SparseArray4i_64 = Array<int64_t, 4, false>;

    /**
     * @brief Sparse one-dimensional array of int64_t values with variable number of components.
     */
    using SparseArrayXi_64 = Array<int64_t, -1, false>;

    /**
     * @brief Sparse one-dimensional array of uint16_t values with one component.
     */
    using SparseArray1ui_16 = Array<uint16_t, 1, false>;

    /**
     * @brief Sparse one-dimensional array of uint16_t values with two components.
     */
    using SparseArray2ui_16 = Array<uint16_t, 2, false>;

    /**
     * @brief Sparse one-dimensional array of uint16_t values with three components.
     */
    using SparseArray3ui_16 = Array<uint16_t, 3, false>;

    /**
     * @brief Sparse one-dimensional array of uint16_t values with four components.
     */
    using SparseArray4ui_16 = Array<uint16_t, 4, false>;

    /**
     * @brief Sparse one-dimensional array of uint16_t values with variable number of components.
     */
    using SparseArrayXui_16 = Array<uint16_t, -1, false>;

    /**
     * @brief Sparse one-dimensional array of uint32_t values with one component.
     */
    using SparseArray1ui = Array<uint32_t, 1, false>;

    /**
     * @brief Sparse one-dimensional array of uint32_t values with two components.
     */
    using SparseArray2ui = Array<uint32_t, 2, false>;

    /**
     * @brief Sparse one-dimensional array of uint32_t values with three components.
     */
    using SparseArray3ui = Array<uint32_t, 3, false>;

    /**
     * @brief Sparse one-dimensional array of uint32_t values with four components.
     */
    using SparseArray4ui = Array<uint32_t, 4, false>;

    /**
     * @brief Sparse one-dimensional array of uint32_t values with variable number of components.
     */
    using SparseArrayXui = Array<uint32_t, -1, false>;

    /**
     * @brief Sparse one-dimensional array of uint64_t values with one component.
     */
    using SparseArray1ui_64 = Array<uint64_t, 1, false>;

    /**
     * @brief Sparse one-dimensional array of uint64_t values with two components.
     */
    using SparseArray2ui_64 = Array<uint64_t, 2, false>;

    /**
     * @brief Sparse one-dimensional array of uint64_t values with three components.
     */
    using SparseArray3ui_64 = Array<uint64_t, 3, false>;

    /**
     * @brief Sparse one-dimensional array of uint64_t values with four components.
     */
    using SparseArray4ui_64 = Array<uint64_t, 4, false>;

    /**
     * @brief Sparse one-dimensional array of uint64_t values with variable number of components.
     */
    using SparseArrayXui_64 = Array<uint64_t, -1, false>;
}
