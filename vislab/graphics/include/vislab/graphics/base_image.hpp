#pragma once

#include "base_image_fwd.hpp"

#include "iimage.hpp"

namespace vislab
{
    /**
     * @brief Abstract base class for two-dimensional images with custom number of color components per pixel.
     * @tparam TChannels Number of color components per pixel.
     * @tparam TD Dummy parameter needed to make this type reflection compatible.
     */
    template <int64_t TChannels>
    class BaseImage : public IImage
    {
        VISLAB_INTERFACE(BaseImage, IImage)

    public:
        /**
         * @brief Gets the number of color channels.
         */
        static const int64_t Channels = TChannels;

        /**
         * @brief Gets the number of channels stored in the image.
         * @return Number of color channels.
         */
        Eigen::Index getNumChannels() const override { return Channels; }

        /**
         * @brief Gets the value at a certain pixel in double precision.
         * @param linearIndex Linear array index.
         * @return Color in double precision.
         */
        virtual Eigen::Vector<double, Channels> getValueDouble(const Eigen::Index& linearIndex) const = 0;

        /**
         * @brief Sets the value at a certain pixel in double precision.
         * @param linearIndex Linear array index.
         * @param value Color in double precision.
         */
        virtual void setValueDouble(const Eigen::Index& linearIndex, const Eigen::Vector<double, Channels>& value) = 0;
    };
}
