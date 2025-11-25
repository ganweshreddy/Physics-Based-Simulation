#include <vislab/graphics/const_texture.hpp>

namespace vislab
{
    ConstTexture::ConstTexture()
        : color(Spectrum(1, 1, 1))
    {
    }

    ConstTexture::ConstTexture(const Spectrum& spectrum)
        : color(spectrum)
    {
    }

    Spectrum ConstTexture::evaluate(const SurfaceInteraction& si) const
    {
        return color;
    }
}
