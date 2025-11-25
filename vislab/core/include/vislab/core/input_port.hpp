#pragma once

#include "iinput_port.hpp"

namespace vislab
{
    /**
     * @brief Strongly typed input port.
     * @tparam TData Data type that is accepted by the port.
     */
    template <typename TData>
    class InputPort : public IInputPort
    {
        VISLAB_OBJECT(InputPort<TData>, IInputPort)

    public:
        /**
         * @brief Gets the data that is stored on the port.
         * @return Shared pointer that might be nullptr if there is no data set.
         */
        std::shared_ptr<const TData> getData() const { return mData; }

        /**
         * @brief Sets the data that is stored on the port.
         * @param data Data to be set on the port.
         */
        void setData(std::shared_ptr<const TData> data) { this->mData = data; }

    private:
        /**
         * @brief Gets the data that is stored on the port.
         * @return Shared pointer that might be nullptr if there is no data set.
         */
        virtual std::shared_ptr<const Data> getDataImpl() const override { return mData; }

    private:
        /**
         * @brief Data stored on the input port.
         */
        std::shared_ptr<const TData> mData;
    };
}
