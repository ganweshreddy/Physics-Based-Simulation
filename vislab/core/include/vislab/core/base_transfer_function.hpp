#pragma once

#include "base_transfer_function_fwd.hpp"
#include "itransfer_function.hpp"

namespace vislab
{
    /**
     * @brief Transfer function with fixed number of output components.
     */
    template <int64_t TComponents>
    class BaseTransferFunction : public ITransferFunction
    {
        VISLAB_INTERFACE(BaseTransferFunction, ITransferFunction)

    public:
        /**
         * @brief Number of components that the transfer function maps to.
         */
        static const int64_t Components = TComponents;

        /**
         * @brief Value type that is mapped to.
         */
        using Value = Eigen::Vector<double, TComponents>;

        /**
         * @brief Constructor. Receives global minimal and maximal transfer function bounds.
         * @param minValue Minimum bound for the transfer function.
         * @param maxValue Maximum bound for the transfer function.
         */
        BaseTransferFunction(const double& minValue = 0, const double& maxValue = 1)
            : ITransferFunction(minValue, maxValue)
        {
        }

        /**
         * @brief Applies the transfer function to a single value.
         * @param value Value from the domain of the transfer function.
         * @return Mapped value.
         */
        virtual Value map(const double& value) const = 0;
    };
}
