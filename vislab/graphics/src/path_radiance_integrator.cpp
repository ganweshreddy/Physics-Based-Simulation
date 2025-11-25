#include <vislab/graphics/path_radiance_integrator.hpp>

#include <vislab/graphics/bsdf.hpp>
#include <vislab/graphics/camera.hpp>
#include <vislab/graphics/image.hpp>
#include <vislab/graphics/light.hpp>
#include <vislab/graphics/sampler.hpp>
#include <vislab/graphics/scene.hpp>

namespace vislab
{
    PathRadianceIntegrator::PathRadianceIntegrator()
        : maxDepth(1)
        , rouletteDepth(1)
    {
    }

    std::pair<Spectrum, bool> PathRadianceIntegrator::sample(const Scene* scene, Sampler* sampler, const RayDifferential3d& ray_) const
    {
        RayDifferential3d ray = ray_;

        // Tracks radiance scaling due to index of refraction changes
        double eta(1.);

        // MIS weight for intersected emitters (set by prev. iteration)
        double emission_weight(1.);

        Spectrum throughput = Spectrum::Ones();
        Spectrum result     = Spectrum::Zero();

        // ---------------------- First intersection ----------------------

        SurfaceInteraction si = scene->intersect(ray);
        bool valid_ray        = si.isValid();
        const Light* light    = si.light(scene);
        bool active           = true;

        for (int64_t depth = 1;; ++depth)
        {
            // --------- hit a light source? ----------
            if (light != nullptr)
                result += emission_weight * throughput.cwiseProduct(light->evaluate(si));

            active &= si.isValid();

            // Early termination with roulette sampling
            if (depth > rouletteDepth)
            {
                double q = std::min(throughput.maxCoeff() * eta * eta, .95);
                active &= sampler->next_1d() < q;
                throughput *= 1. / q;
            }

            // Terminate recursion?
            if (depth > maxDepth || !active)
                break;

            // --------- light sample ----------
            const BSDF* bsdf = si.bsdf(ray);
            bool active_e    = active && bsdf->hasFlag(EBSDFFlag::Smooth);

            if (active_e)
            {
                auto [ds, emitter_val] = scene->sampleLightDirection(
                    si, sampler->next_2d(), true);
                active_e &= ds.pdf != 0.;

                // Query the BSDF for that emitter-sampled direction
                Eigen::Vector3d wo = si.toLocal(ds.direction);
                Spectrum bsdf_val  = bsdf->evaluate(si, wo);

                // Determine density of sampling that same direction using BSDF sampling
                double bsdf_pdf = bsdf->pdf(si, wo);

                double mis = ds.delta ? 1. : misWeight(ds.pdf, bsdf_pdf);
                result += mis * throughput.cwiseProduct(bsdf_val).cwiseProduct(emitter_val);
            }

            // --------- bsdf sample ----------

            // Sample BSDF * cos(theta)
            auto [bs, bsdf_val] = bsdf->sample(si, sampler->next_2d());

            throughput = throughput.cwiseProduct(bsdf_val);
            active &= throughput.any();
            if (!active)
                break;

            eta *= bs.eta;

            // Intersect the BSDF ray against the scene geometry
            ray                        = si.spawnRay(si.toWorld(bs.wo));
            SurfaceInteraction si_bsdf = scene->intersect(ray);

            // Determine probability of having sampled that same direction using emitter sampling.
            light = si_bsdf.light(scene);
            DirectionSample ds(scene, si_bsdf, si);

            if (light != nullptr)
            {
                double emitter_pdf = light != nullptr && !bs.hasFlag(EBSDFFlag::Delta) ? scene->pdfLightDirection(si, ds) : 0.;
                emission_weight    = misWeight(bs.pdf, emitter_pdf);
            }

            si = std::move(si_bsdf);
        }

        return { result, valid_ray };
    }
}
