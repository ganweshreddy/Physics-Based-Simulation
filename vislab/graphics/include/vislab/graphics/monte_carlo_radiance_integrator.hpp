#pragma once

#include "radiance_integrator.hpp"
#include "spectrum.hpp"

#include <memory>

namespace vislab
{
    class Sampler;
    class RayDifferential3d;

    /**
     * @brief Base class for Monte Carlo based integrators.
     */
    class MonteCarloRadianceIntegrator : public RadianceIntegrator
    {
        VISLAB_INTERFACE(MonteCarloRadianceIntegrator, RadianceIntegrator)

    public:
        /**
         * @brief Constructor.
         */
        MonteCarloRadianceIntegrator();

        /**
         * @brief Renders a scene onto a given sensor.
         * @param scene Scene to render.
         * @param sensor Sensor to render scene for.
         * @param image Image to render into.
         * @return True if successful.
         */
        bool render(Scene* scene, Image3d* image) override;

        /**
         * @brief Gets the number of samples per pixel.
         * @return Samples per pixel.
         */
        const size_t& getSamplesPerPixel() const;

        /**
         * @brief Sets the number of samples per pixel.
         * @param spp Samples per pixel.
         */
        void setSamplesPerPixel(const size_t& spp);

        /**
         * @brief Sampler that is cloned for each rendering worker thread.
         */
        std::shared_ptr<Sampler> sampler;

    protected:
        /**
         * @brief Renders a Monte Carlo sample.
         * @param scene Scene to render.
         * @param sampler Sampler used to generate random numbers.
         * @param ray Ray to trace.
         * @return Radiance and true if the sample was valid.
         */
        virtual std::pair<Spectrum, bool> sample(const Scene* scene, Sampler* sampler, const RayDifferential3d& ray) const = 0;

        /**
         * @brief Computes a multiple importance sampling weight using the balance heuristic.
         * @param pdf_a First PDF.
         * @param pdf_b Second PDF.
         * @return MIS weight.
         */
        double misWeight(double pdf_a, double pdf_b) const;

        /**
         * @brief Number of samples per pixel.
         */
        size_t mSamplesPerPixel;
    };
}
