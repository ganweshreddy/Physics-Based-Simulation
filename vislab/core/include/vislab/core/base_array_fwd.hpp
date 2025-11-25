#pragma once

#include <stdint.h>

namespace vislab
{
    /**
     * @brief Abstract base class for data arrays with a fixed number of components.
     * @tparam TDimensions Number of components per element. Use -1 for dynamic number of components.
     */
    template <int64_t TDimensions>
    class BaseArray;

    /**
     * @brief Abstract base array with one component per element.
     */
    using IArray1 = BaseArray<1>;

    /**
     * @brief Abstract base array with two component per element.
     */
    using IArray2 = BaseArray<2>;

    /**
     * @brief Abstract base array with three component per element.
     */
    using IArray3 = BaseArray<3>;

    /**
     * @brief Abstract base array with four component per element.
     */
    using IArray4 = BaseArray<4>;

    /**
     * @brief Abstract base array with variable number of components per element.
     */
    using IArrayX = BaseArray<-1>;
}
