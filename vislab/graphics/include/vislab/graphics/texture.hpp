#pragma once

#include "spectrum.hpp"

#include <vislab/core/data.hpp>

namespace vislab
{
    class SurfaceInteraction;

    /**
     * @brief Base class for textures.
     */
    class Texture : public Data
    {
        VISLAB_INTERFACE(Texture, Data)

    public:
        /**
         * @brief Evaluates the texture at a certain surface interaction location.
         * @param si Surface interaction to evaluate color at.
         * @return Sampled color.
         */
        virtual Spectrum evaluate(const SurfaceInteraction& si) const = 0;

        /**
         * @brief Checks if this texture is spatially-varying.
         * @return True, if spatially-varying.
         */
        virtual bool isSpatiallyVarying() const;
    };
}
