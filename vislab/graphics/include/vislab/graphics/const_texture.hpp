#pragma once

#include "texture.hpp"

namespace vislab
{
    /**
     * @brief Spatially-constant texture with a fixed color.
     */
    class ConstTexture : public Texture
    {
        VISLAB_OBJECT(ConstTexture, Texture)

    public:
        /**
         * @brief Constructor with a white color.
         */
        ConstTexture();

        /**
         * @brief Constructor with an initial color.
         * @param spectrum Initial color.
         */
        ConstTexture(const Spectrum& spectrum);

        /**
         * @brief Evaluates the texture at a certain surface interaction location.
         * @param si Surface interaction to evaluate color at.
         * @return Sampled color.
         */
        Spectrum evaluate(const SurfaceInteraction& si) const override;

        /**
         * @brief Constant color of this texture.
         */
        Spectrum color;
    };
}
