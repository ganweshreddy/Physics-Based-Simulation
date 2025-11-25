#pragma once

#include "monte_carlo_radiance_integrator.hpp"

namespace vislab
{
    /**
     * @brief Simple direct radiance integrator.
     */
    class DirectRadianceIntegrator : public MonteCarloRadianceIntegrator
    {
        VISLAB_OBJECT(DirectRadianceIntegrator, MonteCarloRadianceIntegrator)

    public:
        /**
         * @brief Constructor.
         */
        DirectRadianceIntegrator();

        /**
         * @brief Renders a Monte Carlo sample.
         * @param scene Scene to render.
         * @param sampler Sampler used to generate random numbers.
         * @param ray Ray to trace.
         * @return Radiance and true if the sample was valid.
         */
        std::pair<Spectrum, bool> sample(const Scene* scene, Sampler* sampler, const RayDifferential3d& ray) const override;

        /**
         * @brief Number of light samples.
         */
        int64_t lightSamples;

        /**
         * @brief Number of bsdf samples.
         */
        int64_t bsdfSamples;
    };
}
