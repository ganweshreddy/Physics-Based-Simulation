#pragma once

#include "data.hpp"

namespace vislab
{
    /**
     * @brief Interface for transfer functions.
     */
    class ITransferFunction : public Data
    {
        VISLAB_INTERFACE(ITransferFunction, Data)

    public:
        /**
         * @brief Constructor. Receives global minimal and maximal transfer function bounds.
         * @param minValue Minimum bound for the transfer function.
         * @param maxValue Maximum bound for the transfer function.
         */
        ITransferFunction(const double& minValue = 0, const double& maxValue = 1);

        /**
         * @brief Global minimal bounds of the transfer function.
         */
        double minValue;

        /**
         * @brief Global maximal bounds of the transfer function.
         */
        double maxValue;
    };
}
