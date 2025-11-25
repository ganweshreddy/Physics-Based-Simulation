#pragma once

#include <vislab/core/data.hpp>
#include <vislab/core/types.hpp>

namespace vislab
{
    /**
     * @brief Interface for image data structures in any dimension.
     */
    class IImage : public Data
    {
        VISLAB_INTERFACE(IImage, Data)

    public:
        /**
         * @brief Gets the total number of pixels.
         * @return Number of pixels.
         */
        virtual Eigen::Index getNumPixels() const = 0;

        /**
         * @brief Gets the number of channels stored in the image.
         * @return Number of color channels.
         */
        virtual Eigen::Index getNumChannels() const = 0;

        /**
         * @brief Gets the resolution of the image.
         * @return Image resolution.
         */
        virtual const Eigen::Vector2i& getResolution() const = 0;

        /**
         * @brief Sets the resolution of the image.
         * @param resolution Image resolution to set.
         */
        virtual void setResolution(const Eigen::Vector2i& resolution) = 0;

        /**
         * @brief Sets the resolution of the image.
         * @param resx New x resolution.
         * @param resy New y resolution.
         */
        virtual void setResolution(const int resx, const int resy) = 0;

        /**
         * @brief Sets all values to zero.
         */
        virtual void setZero() = 0;
    };
}
