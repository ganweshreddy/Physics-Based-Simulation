#include <vislab/graphics/monte_carlo_radiance_integrator.hpp>

#include <vislab/graphics/camera.hpp>
#include <vislab/graphics/image.hpp>
#include <vislab/graphics/independent_sampler.hpp>
#include <vislab/graphics/ray.hpp>
#include <vislab/graphics/sampler.hpp>
#include <vislab/graphics/scene.hpp>
#include <vislab/graphics/spectrum.hpp>

namespace vislab
{
    MonteCarloRadianceIntegrator::MonteCarloRadianceIntegrator()
        : mSamplesPerPixel(1)
        , sampler(std::make_shared<IndependentSampler>())
    {
    }

    bool MonteCarloRadianceIntegrator::render(Scene* scene, Image3d* image)
    {
        // get the camera to render from
        auto camera = scene->camera;

        // for each pixel in the image
        Eigen::Vector2i resolution = image->getResolution();

#ifndef _DEBUG
#pragma omp parallel for
#endif
        for (int64_t i = 0; i < (int64_t)resolution.x() * resolution.y(); ++i)
        {
            // make a copy of the random sampler, such that each pixel has its own
            std::shared_ptr<IndependentSampler> s = std::dynamic_pointer_cast<IndependentSampler>(sampler);
            std::shared_ptr<Sampler> local_sampler(s->clone());
            local_sampler->seed(i);

            // determine the pixel index
            Eigen::Vector2d pixelIndex = Eigen::Vector2d(
                i % resolution.x(),
                i / resolution.x());

            // iterate to obtain the requested number of samples per pixel
            size_t numSamples = 0;
            Spectrum radiance = Spectrum::Zero();
            for (size_t spp = 0; spp < mSamplesPerPixel; ++spp)
            {
                // generate a random sample inside the pixel
                Eigen::Vector2d s = (pixelIndex + local_sampler->next_2d()).cwiseQuotient(resolution.cast<double>());

                // generate the ray
                RayDifferential3d ray;
                Spectrum spec_weight;
                std::tie(ray, spec_weight) = camera->sampleRayDifferential(s);
                ray.medium                 = camera->medium.get();

                // let the integrator compute the radiance
                Spectrum spectrum;
                bool isValid;
                std::tie(spectrum, isValid) = this->sample(scene, local_sampler.get(), ray);

                // if valid sample
                if (isValid)
                {
                    // add the color to the buffer
                    radiance += spectrum;
                    numSamples++;
                }
            }
            // divide by the total number of samples
            if (numSamples > 0)
            {
                radiance /= numSamples;
                // apply gamma correction
                radiance.x() = std::min(std::max(0., std::pow(radiance.x(), 1. / 2.2)), 1.);
                radiance.y() = std::min(std::max(0., std::pow(radiance.y(), 1. / 2.2)), 1.);
                radiance.z() = std::min(std::max(0., std::pow(radiance.z(), 1. / 2.2)), 1.);
                image->setValue(i, radiance);
            }
        }
        return true;
    }

    const size_t& MonteCarloRadianceIntegrator::getSamplesPerPixel() const { return mSamplesPerPixel; }
    void MonteCarloRadianceIntegrator::setSamplesPerPixel(const size_t& spp) { mSamplesPerPixel = spp; }

    double MonteCarloRadianceIntegrator::misWeight(double pdf_a, double pdf_b) const
    {
        pdf_a *= pdf_a;
        pdf_b *= pdf_b;
        double w = pdf_a / (pdf_a + pdf_b);
        return w == w ? w : 0;
    }
}
