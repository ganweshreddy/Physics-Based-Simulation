#include <vislab/graphics/dielectric_bsdf.hpp>

#include <vislab/graphics/const_texture.hpp>
#include <vislab/graphics/fresnel.hpp>
#include <vislab/graphics/texture.hpp>
#include <vislab/graphics/warp.hpp>

namespace vislab
{
    DielectricBSDF::DielectricBSDF()
        : BSDF((uint32_t)EBSDFFlag::DeltaReflection | (uint32_t)EBSDFFlag::DeltaTransmission | (uint32_t)EBSDFFlag::FrontSide | (uint32_t)EBSDFFlag::BackSide | (uint32_t)EBSDFFlag::NonSymmetric)
        , interior_ior(1.5046)
        , exterior_ior(1.000277)
        , eta(1.5046 / 1.000277)
        , specularReflectance(std::make_shared<ConstTexture>(Spectrum(1.0, 1.0, 1.0)))
        , specularTransmittance(std::make_shared<ConstTexture>(Spectrum(1.0, 1.0, 1.0)))
    {
    }

    std::pair<BSDFSample, Spectrum> DielectricBSDF::sample(const SurfaceInteraction& si, const Eigen::Vector2d& sample) const
    {
        // Evaluate the Fresnel equations for unpolarized illumination
        double cos_theta_i = Frame::cos_theta(si.wi);

        auto [r_i, cos_theta_t, eta_it, eta_ti] = fresnel(cos_theta_i, eta);
        double t_i                              = 1.f - r_i;

        bool selected_r = sample.x() <= r_i;
        bool selected_t = !selected_r;

        BSDFSample bs(Eigen::Vector3d(0, 0, 0));
        bs.pdf         = selected_r ? r_i : t_i;
        bs.sampledType = selected_r ? uint32_t(EBSDFFlag::DeltaReflection) : uint32_t(EBSDFFlag::DeltaTransmission);
        bs.wo          = selected_r ? reflect(si.wi) : refract(si.wi, cos_theta_t, eta_ti);
        bs.eta         = selected_r ? 1. : eta_it;

        Spectrum reflectance = Spectrum::Ones(), transmittance = Spectrum::Ones();
        if (specularReflectance && selected_r)
            reflectance = specularReflectance->evaluate(si);
        if (specularTransmittance && selected_t)
            transmittance = specularTransmittance->evaluate(si);

        Spectrum weight = Spectrum::Ones();
        if (selected_r)
            weight = weight.cwiseProduct(reflectance);
        if (selected_t)
        {
            weight = weight.cwiseProduct(transmittance);
            // For transmission, radiance must be scaled to account for the solid angle compression that occurs when crossing the interface.
            double factor = eta_ti;
            weight *= (factor * factor);
        }
        return { bs, weight };
    }

    Spectrum DielectricBSDF::evaluate(const SurfaceInteraction& si, const Eigen::Vector3d& wo) const
    {
        return Spectrum::Zero();
    }

    double DielectricBSDF::pdf(const SurfaceInteraction& si, const Eigen::Vector3d& wo) const
    {
        return 0.;
    }
}
