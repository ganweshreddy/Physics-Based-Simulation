#include <vislab/graphics/colormap_texture.hpp>

#include <vislab/graphics/surface_interaction.hpp>

#include <vislab/field/regular_field.hpp>

namespace vislab
{
    ColormapTexture::ColormapTexture()
    {
        transferFunction.minValue = 0;
        transferFunction.maxValue = 1;
        transferFunction.values.insert(std::make_pair(0, Eigen::Vector4d::Zero()));
        transferFunction.values.insert(std::make_pair(1, Eigen::Vector4d::Ones()));
    }

    Spectrum ColormapTexture::evaluate(const SurfaceInteraction& si) const
    {
        Eigen::AlignedBox2d domain = scalarField->getDomain();
        Eigen::Vector2d pos        = domain.min() + (domain.max() - domain.min()).cwiseProduct(si.uv);
        double scalarValue         = scalarField->sample(pos).x();
        Eigen::Vector4d color      = transferFunction.map(scalarValue);
        return color.xyz();
    }
}
