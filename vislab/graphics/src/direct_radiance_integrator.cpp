#include <vislab/graphics/direct_radiance_integrator.hpp>

#include <vislab/graphics/bsdf.hpp>
#include <vislab/graphics/camera.hpp>
#include <vislab/graphics/image.hpp>
#include <vislab/graphics/light.hpp>
#include <vislab/graphics/sampler.hpp>
#include <vislab/graphics/scene.hpp>

namespace vislab
{
    DirectRadianceIntegrator::DirectRadianceIntegrator()
        : lightSamples(1)
        , bsdfSamples(1)
    {
    }

    std::pair<Spectrum, bool> DirectRadianceIntegrator::sample(const Scene* scene, Sampler* sampler, const RayDifferential3d& ray) const
    {
        // shoot ray into the scene
        auto si = scene->intersect(ray);

        Spectrum result = Spectrum::Zero();

        // --------- visible emitters ----------

        // if hit a light, evaluate it
        const Light* vis_light = si.light(scene);
        if (vis_light)
        {
            result += vis_light->evaluate(si);
        }
        if (!si.isValid())
        {
            return std::make_pair(result, false);
        }

        double sum         = lightSamples + bsdfSamples;
        double bsdfWeight  = 1. / bsdfSamples;
        double lightWeight = 1. / lightSamples;
        double bsdfFrac    = bsdfSamples / sum;
        double lightFrac   = lightSamples / sum;

        // --------- light samples ----------
        const BSDF* bsdf   = si.bsdf(ray);
        bool sampleEmitter = bsdf->hasFlag(EBSDFFlag::Smooth);

        if (sampleEmitter)
        {
            for (int i = 0; i < lightSamples; ++i)
            {
                // sample the light source
                DirectionSample ds;
                Spectrum emitter_val;
                std::tie(ds, emitter_val) = scene->sampleLightDirection(si, sampler->next_2d(), true);
                if (ds.pdf == 0)
                    continue;

                // Query the BSDF for that emitter-sampled direction
                Eigen::Vector3d wo = si.toLocal(ds.direction);

                // Determine BSDF value and probability of having sampled that same direction using BSDF sampling.
                auto [bsdf_val, bsdf_pdf] = bsdf->evaluate_pdf(si, wo);

                double mis = ds.delta ? 1. : misWeight(ds.pdf * lightFrac, bsdf_pdf * bsdfFrac) * lightWeight;
                result += mis * bsdf_val.cwiseProduct(emitter_val);
            }
        }

        // --------- bsdf samples ----------
        for (int i = 0; i < bsdfSamples; ++i)
        {
            auto [bs, bsdf_val] = bsdf->sample(si, sampler->next_2d());
            // bsdf_val            = si.to_world_mueller(bsdf_val, -bs.wo, si.wi);

            bool active_b = bsdf_val.any();

            // Trace the ray in the sampled direction and intersect against the scene
            SurfaceInteraction si_bsdf =
                scene->intersect(si.spawnRay(si.toWorld(bs.wo)));

            // Retain only rays that hit an emitter
            const Light* light = si_bsdf.light(scene);
            active_b &= light != nullptr;

            if (active_b)
            {
                Spectrum light_val = light->evaluate(si_bsdf);
                bool delta         = bs.hasFlag(EBSDFFlag::Delta);

                // Determine probability of having sampled that same direction using Emitter sampling.
                DirectionSample ds(scene, si_bsdf, si);

                double light_pdf = delta ? 0.f : scene->pdfLightDirection(si, ds);

                double mis = misWeight(bs.pdf * bsdfFrac, light_pdf * lightFrac) * bsdfWeight;
                result += mis * bsdf_val.cwiseProduct(light_val);
            }
        }

        return std::make_pair(result, true);
    }
}
