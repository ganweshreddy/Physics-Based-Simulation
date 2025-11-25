#pragma once

#include "data.hpp"

namespace vislab
{
    /**
     * @brief Basic interface for a parameter.
     */
    class IParameter : public Data
    {
        VISLAB_INTERFACE(IParameter, Data)
    };
}
