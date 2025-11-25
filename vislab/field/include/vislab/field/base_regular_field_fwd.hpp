#pragma once

#include "regular_grid_fwd.hpp"

namespace vislab
{
    /**
     * @brief Base class for fields defined on regular grids.
     * @tparam TBaseType Field type to discretize onto a regular grid.
     * @tparam TGridType Grid type to discretize on.
     */
    template <class TBaseType, typename TGridType>
    class BaseRegularField;

    /**
     * @brief Base class for 2D steady scalar fields.
     */
    class ISteadyScalarField2d;

    /**
     * @brief Base class for 3D steady scalar fields.
     */
    class ISteadyScalarField3d;

    /**
     * @brief Base class for 2D unsteady scalar fields.
     */
    class IUnsteadyScalarField2d;

    /**
     * @brief Base class for 3D unsteady scalar fields.
     */
    class IUnsteadyScalarField3d;

    /**
     * @brief Base class for 2D steady vector fields.
     */
    class ISteadyVectorField2d;

    /**
     * @brief Base class for 3D steady vector fields.
     */
    class ISteadyVectorField3d;

    /**
     * @brief Base class for 2D unsteady vector fields.
     */
    class IUnsteadyVectorField2d;

    /**
     * @brief Base class for 3D unsteady vector fields.
     */
    class IUnsteadyVectorField3d;

    /**
     * @brief Two-dimensional steady scalar field on a regular grid.
     */
    using IRegularSteadyScalarField2 = BaseRegularField<ISteadyScalarField2d, RegularGrid2d>;

    /**
     * @brief Three-dimensional steady scalar field on a regular grid.
     */
    using IRegularSteadyScalarField3 = BaseRegularField<ISteadyScalarField3d, RegularGrid3d>;

    /**
     * @brief Two-dimensional unsteady scalar field on a regular grid.
     */
    using IRegularUnsteadyScalarField2 = BaseRegularField<IUnsteadyScalarField2d, RegularGrid3d>;

    /**
     * @brief Three-dimensional unsteady scalar field on a regular grid.
     */
    using IRegularUnsteadyScalarField3 = BaseRegularField<IUnsteadyScalarField3d, RegularGrid4d>;

    /**
     * @brief Two-dimensional steady vector field on a regular grid.
     */
    using IRegularSteadyVectorField2 = BaseRegularField<ISteadyVectorField2d, RegularGrid2d>;

    /**
     * @brief Three-dimensional steady vector field on a regular grid.
     */
    using IRegularSteadyVectorField3 = BaseRegularField<ISteadyVectorField3d, RegularGrid3d>;

    /**
     * @brief Two-dimensional unsteady vector field on a regular grid.
     */
    using IRegularUnsteadyVectorField2 = BaseRegularField<IUnsteadyVectorField2d, RegularGrid3d>;

    /**
     * @brief Three-dimensional unsteady vector field on a regular grid.
     */
    using IIRegularUnsteadyVectorField3 = BaseRegularField<IUnsteadyVectorField3d, RegularGrid4d>;
}
