#pragma once

#include "data.hpp"

#include <memory>

namespace vislab
{
    /**
     * @brief Interface for input ports.
     */
    class IInputPort : public Object
    {
        VISLAB_INTERFACE(IInputPort, Object)

    public:
        /**
         * @brief Constructor.
         * @param isRequired Optional flag that sets whether the input port is required for execution.
         */
        IInputPort(bool isRequired = true)
            : required(isRequired)
        {
        }

        /**
         * @brief Gets the data that is stored on the port.
         * @return Shared pointer that might be nullptr if there is no data set.
         */
        inline std::shared_ptr<const Data> getData() const { return getDataImpl(); }

        /**
         * @brief Flag that determines whether the input port is required for executing the algorithm.
         */
        bool required;

    private:
        /**
         * @brief Gets the data that is stored on the port.
         * @return Shared pointer that might be nullptr if there is no data set.
         */
        virtual std::shared_ptr<const Data> getDataImpl() const = 0;
    };
}
