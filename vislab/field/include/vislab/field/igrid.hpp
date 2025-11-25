#pragma once

#include <vislab/core/data.hpp>
#include <vislab/core/types.hpp>

namespace vislab
{
    /**
     * @brief Basic interface for grid data structures in any dimension.
     */
    class IGrid : public Data
    {
        VISLAB_INTERFACE(IGrid)

    public:
        /**
         * @brief Gets the total number of grid points.
         * @return Number of grid points.
         */
        virtual Eigen::Index getNumGridPoints() const = 0;

        /**
         * @brief Gets the number of cells.
         * @return Number of cells.
         */
        virtual Eigen::Index getNumCells() const = 0;

        /**
         * @brief Gets the linear grid point indices of a given cell.
         * @param cellIndex Index of cell to get.
         * @return Grid point indices.
         */
        virtual Eigen::VectorX<Eigen::Index> getCell(const Eigen::Index& cellIndex) const = 0;
    };
}
