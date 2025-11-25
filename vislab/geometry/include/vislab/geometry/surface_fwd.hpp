#pragma once

#include <vislab/core/array_fwd.hpp>

namespace vislab
{
    /**
     * @brief Surface geometry with specific internal type.
     * @tparam TArrayType Internal array type that stores the vertex positions.
     */
    template <typename TArrayType>
    class Surface;

    /**
     * @brief One-dimensional surface geometry with vertex positions in float precision.
     */
    using Surface1f = Surface<Array1f>;

    /**
     * @brief Two-dimensional surface geometry with vertex positions in float precision.
     */
    using Surface2f = Surface<Array2f>;

    /**
     * @brief Three-dimensional surface geometry with vertex positions in float precision.
     */
    using Surface3f = Surface<Array3f>;

    /**
     * @brief Four-dimensional surface geometry with vertex positions in float precision.
     */
    using Surface4f = Surface<Array4f>;

    /**
     * @brief One-dimensional surface geometry with vertex positions in double precision.
     */
    using Surface1d = Surface<Array1d>;

    /**
     * @brief Two-dimensional surface geometry with vertex positions in double precision.
     */
    using Surface2d = Surface<Array2d>;

    /**
     * @brief Three-dimensional surface geometry with vertex positions in double precision.
     */
    using Surface3d = Surface<Array3d>;

    /**
     * @brief Four-dimensional surface geometry with vertex positions in double precision.
     */
    using Surface4d = Surface<Array4d>;
}
