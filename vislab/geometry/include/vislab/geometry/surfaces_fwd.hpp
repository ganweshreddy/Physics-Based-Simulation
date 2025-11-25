#pragma once

#include <vislab/core/array_fwd.hpp>

namespace vislab
{
    /**
     * @brief Set of surface geometries with specific internal type.
     * @tparam TArrayType Internal array type that stores the vertex positions.
     */
    template <typename TArrayType>
    class Surfaces;

    /**
     * @brief One-dimensional sets of surface geometries with vertex positions in float precision.
     */
    using Surfaces1f = Surfaces<Array1f>;

    /**
     * @brief Two-dimensional sets of surface geometries with vertex positions in float precision.
     */
    using Surfaces2f = Surfaces<Array2f>;

    /**
     * @brief Three-dimensional sets of surface geometries with vertex positions in float precision.
     */
    using Surfaces3f = Surfaces<Array3f>;

    /**
     * @brief Four-dimensional sets of surface geometries with vertex positions in float precision.
     */
    using Surfaces4f = Surfaces<Array4f>;

    /**
     * @brief One-dimensional sets of surface geometries with vertex positions in double precision.
     */
    using Surfaces1d = Surfaces<Array1d>;

    /**
     * @brief Two-dimensional sets of surface geometries with vertex positions in double precision.
     */
    using Surfaces2d = Surfaces<Array2d>;

    /**
     * @brief Three-dimensional sets of surface geometries with vertex positions in double precision.
     */
    using Surfaces3d = Surfaces<Array3d>;

    /**
     * @brief Four-dimensional sets of surface geometries with vertex positions in double precision.
     */
    using Surfaces4d = Surfaces<Array4d>;
}
