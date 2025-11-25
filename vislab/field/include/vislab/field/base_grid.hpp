#pragma once

#include "base_grid_fwd.hpp"

#include "igrid.hpp"

namespace vislab
{
    /**
     * @brief Interface for grid data structures with a certain dimension.
     * @tparam TDimensions Total number of dimensions.
     */
    template <int64_t TDimensions>
    class BaseGrid : public IGrid
    {
        VISLAB_INTERFACE(BaseGrid, IGrid)

    public:
        /**
         * @brief Gets the number of dimensions.
         */
        static const int64_t Dimensions = TDimensions;

        /**
         * @brief Gets the physical coordinate of a grid point in double precision, identified by its linear index.
         * @param linearIndex Linear index of grid point.
         * @return Physical domain coordinates of grid point.
         */
        virtual Eigen::Vector<double, TDimensions> getCoordAtDouble(const Eigen::Index& linearIndex) const = 0;
    };
}
