#pragma once

#include "base_image.hpp"

#include <vislab/core/algorithm.hpp>
#include <vislab/core/path_parameter.hpp>

namespace vislab
{
    /**
     * @brief Bitmap writer for RGB images.
     */
    class BmpWriter : public Algorithm
    {
        VISLAB_ALGORITHM(BmpWriter, Algorithm)

    public:
        /**
         * @brief Constructor.
         */
        BmpWriter();

        /**
         * @brief Image data to write to a bitmap file.
         */
        InputPort<IImage3> inputImage;

        /**
         * @brief Path to write to.
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
