#pragma once

#include "monte_carlo_radiance_integrator.hpp"

namespace vislab
{
    /**
     * @brief Path tracing radiance integrator.
     */
    class PathRadianceIntegrator : public MonteCarloRadianceIntegrator
    {
        VISLAB_OBJECT(PathRadianceIntegrator, MonteCarloRadianceIntegrator)

    public:
        /**
         * @brief Constructor.
         */
        PathRadianceIntegrator();

        /**
         * @brief Renders a Monte Carlo sample.
         * @param scene Scene to render.
         * @param sampler Sampler used to generate random numbers.
         * @param ray Ray to trace.
         * @return Radiance and true if the sample was valid.
         */
        std::pair<Spectrum, bool> sample(const Scene* scene, Sampler* sampler, const RayDifferential3d& ray) const override;

        /**
         * @brief Maximum recursive path tracing depth.
         */
        int64_t maxDepth;

        /**
         * @brief Roulette termination begins at this depth.
         */
        int64_t rouletteDepth;
    };
}
