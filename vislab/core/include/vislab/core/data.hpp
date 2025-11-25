#pragma once

#include "object.hpp"

namespace vislab
{
    /**
     * @brief Base class for a data object.
     */
    class Data : public Object
    {
        VISLAB_INTERFACE(Data, Object)
    };
}
