#pragma once

#include "data.hpp"

#include <memory>

namespace vislab
{
    /**
     * @brief Interface for output ports.
     */
    class IOutputPort : public Object
    {
        VISLAB_INTERFACE(IOutputPort, Object)

    public:
        /**
         * @brief Gets the data that is stored on the port.
         * @return Shared pointer that might be nullptr if there is no data.
         */
        inline std::shared_ptr<Data> getData() { return getDataImpl(); }

    private:
        /**
         * @brief Gets the data that is stored on the port.
         * @return Shared pointer that might be nullptr if there is no data set.
         */
        virtual std::shared_ptr<Data> getDataImpl() = 0;
    };
}
