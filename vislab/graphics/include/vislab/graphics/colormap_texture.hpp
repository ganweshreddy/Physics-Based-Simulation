#pragma once

#include "texture.hpp"

#include <memory>
#include <vislab/core/event.hpp>
#include <vislab/core/transfer_function.hpp>
#include <vislab/field/regular_field_fwd.hpp>

namespace vislab
{
    /**
     * @brief Texture that is derived via a transfer function from a scalar field.
     */
    class ColormapTexture : public Texture
    {
        VISLAB_OBJECT(ColormapTexture, Texture)

    public:
        /**
         * @brief Constructor a default color map from [0,1] -> black-white
         */
        ColormapTexture();

        /**
         * @brief Evaluates the texture at a certain surface interaction location.
         * @param si Surface interaction to evaluate color at.
         * @return Sampled color.
         */
        Spectrum evaluate(const SurfaceInteraction& si) const override;

        /**
         * @brief Transfer function that maps the scalar value to color.
         */
        TransferFunction4d transferFunction;

        /**
         * @brief Scalar field that is mapped to color.
         */
        std::shared_ptr<RegularSteadyScalarField2f> scalarField;

        /**
         * @brief Event that can be raised when the positions changed.
         */
        TEvent<ColormapTexture, void> scalarFieldChanged;
    };
}
