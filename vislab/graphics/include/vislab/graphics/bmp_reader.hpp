#pragma once

#include "image.hpp"

#include <vislab/core/algorithm.hpp>
#include <vislab/core/path_parameter.hpp>

namespace vislab
{
    /**
     * @brief Bitmap reader for RGB images.
     */
    class BmpReader : public Algorithm
    {
        VISLAB_ALGORITHM(BmpReader, Algorithm)

    public:
        /**
         * @brief Constructor.
         */
        BmpReader();

        /**
         * @brief Image data that was read from file.
         */
        OutputPort<Image3d> outputImage;

        /**
         * @brief Path to read from.
         */
        PathParameter paramPath;

    protected:
        /**
         * @brief Internal computation function
         * @param progress Optional progress info.
         * @return Information about the completion of the computation, including a potential error message.
         */
        UpdateInfo internalUpdate(ProgressInfo& progress) override;
    };
}
