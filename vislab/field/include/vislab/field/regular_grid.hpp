#pragma once

#include "base_grid.hpp"
#include "regular_grid_fwd.hpp"

namespace vislab
{
    /**
     * @brief Basic interface for regular grids in a certain dimension.
     * @tparam TScalar Scalar value used to store data on the grid.
     * @tparam TDimensions Number of dimensions.
     */
    template <typename TScalar, int64_t TDimensions>
    class RegularGrid : public BaseGrid<TDimensions>
    {
        VISLAB_OBJECT(RegularGrid, BaseGrid<TDimensions>)

    public:
        /**
         * @brief Number of dimensions.
         */
        static const int64_t Dimensions = TDimensions;

        /**
         * @brief Scalar type of domain coordinates.
         */
        using Scalar = TScalar;

        /**
         * @brief Grid index coordinate.
         */
        using GridCoord = typename Eigen::Vector<int, TDimensions>;

        /**
         * @brief Axis-aligned bounding box type.
         */
        using BoundingBox = typename Eigen::AlignedBox<TScalar, TDimensions>;

        /**
         * @brief Physical domain coordinate.
         */
        using DomainCoord = typename Eigen::Vector<TScalar, TDimensions>;

        /**
         * @brief Gets the grid resolution.
         * @return Resolution of the grid.
         */
        const GridCoord& getResolution() const { return this->mResolution; }

        /**
         * @brief Sets the grid resolution.
         * @param resolution Resolution to set.
         */
        void setResolution(const GridCoord& resolution) { this->mResolution = resolution; }

        /**
         * @brief Gets the domain bounding box.
         * @return Bounding box of domain.
         */
        const BoundingBox& getDomain() const { return this->mDomain; }

        /**
         * @brief Sets the domain bounding box.
         * @param domain Bounding box to set.
         */
        void setDomain(const BoundingBox& domain) { this->mDomain = domain; }

        /**
         * @brief Gets the total number of grid points.
         * @return Total number of grid points.
         */
        Eigen::Index getNumGridPoints() const override
        {
            return mResolution.prod();
        }

        /**
         * @brief Gets the location of a grid point in the physical domain in double precision.
         * @param linearIndex Linear index of grid point.
         * @return Physical domain coordinate of grid point.
         */
        Eigen::Vector<double, Dimensions> getCoordAtDouble(const Eigen::Index& linearIndex) const override
        {
            return getCoordAt(linearIndex).template cast<double>();
        }

        /**
         * @brief Gets the physical coordinate of a grid point, identified by its linear index.
         * @param linearIndex Linear index of grid point.
         * @return Physical domain coordinate of grid point.
         */
        virtual DomainCoord getCoordAt(const Eigen::Index& linearIndex) const
        {
            return getCoordAt(getGridCoord(linearIndex));
        }

        /**
         * @brief Gets the linear array index based on a grid coordinate index.
         * @param gridCoord Grid coordinate.
         * @return Corresponding linear index.
         */
        Eigen::Index getLinearIndex(const GridCoord& gridCoord) const
        {
            Eigen::Index stride      = 1;
            Eigen::Index linearIndex = gridCoord[0];
            for (int64_t d = 1; d < Dimensions; ++d)
            {
                stride *= this->mResolution[d - 1];
                linearIndex += gridCoord[d] * stride;
            }
            return linearIndex;
        }

        /**
         * @brief Gets the spatial location of a grid vertex.
         * @param gridCoord Grid coordinate.
         * @return Physical domain coordiante of grid point.
         */
        DomainCoord getCoordAt(const GridCoord& gridCoord) const
        {
            DomainCoord s;
            for (int i = 0; i < Dimensions; ++i)
            {
                s[i] = mResolution[i] < 2 ? 0.5 : gridCoord[i] / (this->mResolution[i] - Scalar(1.));
            }
            return this->mDomain.min() + (this->mDomain.max() - this->mDomain.min()).cwiseProduct(s);
        }

        /**
         * @brief Gets the grid coordinate based on the linear array index.
         * @param linearIndex Linear array index.
         * @return Corresponding grid coordinate.
         */
        GridCoord getGridCoord(const Eigen::Index& linearIndex) const
        {
            GridCoord result;
            Eigen::Index stride = 1;
            for (int d = 0; d < Dimensions - 1; ++d)
                stride *= this->mResolution[d];

            Eigen::Index t = linearIndex;
            for (int64_t d = Dimensions - 1; d >= 0; --d)
            {
                result[d] = (int)(t / stride);
                t         = t % stride;
                if (d > 0)
                    stride /= this->mResolution[d - 1];
            }
            return result;
        }

        /**
         * @brief Gets the number of cells.
         * @return Number of cells.
         */
        Eigen::Index getNumCells() const override
        {
            Eigen::Index points = 1;
            for (int i = 0; i < Dimensions; ++i)
                points *= (Eigen::Index)this->mResolution(i) - 1;
            return points;
        }

        /**
         * @brief Gets the linear grid point indices of a given cell.
         * @param cellIndex Index of the cell to get.
         * @return Grid point coordinates of the cell.
         */
        Eigen::VectorX<Eigen::Index> getCell(const Eigen::Index& cellIndex) const override
        {
            throw std::logic_error("Function not yet implemented");
        }

        /**
         * @brief Gets the origin of the domain.
         * @return Origin of the domain.
         */
        DomainCoord getOrigin() const { return mDomain.min(); }

        /**
         * @brief Gets the spacing between adjacent grid points.
         * @return Spacing of the domain.
         */
        DomainCoord getSpacing() const { return (mDomain.max() - mDomain.min()).cwiseQuotient((mResolution - GridCoord::Ones()).template cast<double>()); }

    private:
        /**
         * @brief Domain over which the regular grid is defined.
         */
        BoundingBox mDomain;

        /**
         * @brief Resolution of this regular grid.
         */
        GridCoord mResolution;
    };
}
