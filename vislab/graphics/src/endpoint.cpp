#include <vislab/graphics/endpoint.hpp>

namespace vislab
{
    std::pair<DirectionSample, Spectrum> Endpoint::sampleDirection(const Interaction& it, const Eigen::Vector2d& sample) const
    {
        throw std::logic_error("Not implemented.");
        return { DirectionSample(), Spectrum::Zero() };
    }

    double Endpoint::pdfDirection(const Interaction& it, const DirectionSample& ds) const
    {
        throw std::logic_error("Not implemented.");
        return 0;
    }

    Spectrum Endpoint::evaluate(const SurfaceInteraction& si) const
    {
        throw std::logic_error("Not implemented.");
        return Spectrum::Zero();
    }
}
