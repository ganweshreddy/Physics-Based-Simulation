#pragma once

#include <vislab/core/data.hpp>
#include <vislab/core/types.hpp>

namespace vislab
{
    /**
     * @brief Basic interface for any type of field.
     */
    class IField : public Data
    {
        VISLAB_INTERFACE(IField)

    public:
        /**
         * @brief Gets the number of input dimensions to this field (space + time).
         * @return Number of dimensions.
         */
        virtual Eigen::Index getDimensions() const = 0;

        /**
         * @brief Gets the number of input dimensions to this field (space only).
         * @return Number of spatial dimensions.
         */
        virtual Eigen::Index getSpatialDimensions() const = 0;
    };
}
