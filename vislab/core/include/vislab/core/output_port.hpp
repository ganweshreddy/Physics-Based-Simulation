#pragma once

#include "ioutput_port.hpp"

namespace vislab
{
    /**
     * @brief Strongly typed output port.
     * @tparam TData Data type that is accepted by the port.
     */
    template <typename TData>
    class OutputPort : public IOutputPort
    {
        VISLAB_OBJECT(OutputPort<TData>, IOutputPort)

    public:
        /**
         * @brief Gets the data that is stored on the port.
         * @return Weak pointer that might be nullptr if there is no data set.
         */
        std::shared_ptr<TData> getData() { return mData; }

        /**
         * @brief Sets the data that is stored on the port.
         * @param data Data to be set on the port.
         */
        void setData(std::shared_ptr<TData> data) { this->mData = data; }

    private:
        /**
         * @brief Gets the data that is stored on the port.
         * @return Weak pointer that might be nullptr if there is no data set.
         */
        virtual std::shared_ptr<Data> getDataImpl() override { return mData; }

    private:
        /**
         * @brief Data stored on the input port.
         */
        std::shared_ptr<TData> mData;
    };
}
