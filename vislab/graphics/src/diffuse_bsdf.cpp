#include <vislab/graphics/diffuse_bsdf.hpp>

#include <vislab/graphics/const_texture.hpp>
#include <vislab/graphics/texture.hpp>
#include <vislab/graphics/warp.hpp>

namespace vislab
{
    DiffuseBSDF::DiffuseBSDF()
        : BSDF((uint32_t)EBSDFFlag::DiffuseReflection | (uint32_t)EBSDFFlag::FrontSide)
        , reflectance(std::make_shared<ConstTexture>(Spectrum(0.5, 0.5, 0.5)))
    {
    }

    DiffuseBSDF::DiffuseBSDF(std::shared_ptr<Texture> reflectance)
        : BSDF((uint32_t)EBSDFFlag::DiffuseReflection | (uint32_t)EBSDFFlag::FrontSide)
        , reflectance(reflectance)
    {
    }

    std::pair<BSDFSample, Spectrum> DiffuseBSDF::sample(const SurfaceInteraction& si, const Eigen::Vector2d& sample) const
    {
        double cos_theta_i = Frame::cos_theta(si.wi);
        BSDFSample bs(Eigen::Vector3d::Zero());

        bool active = cos_theta_i > 0.;
        if (!active)
        {
            return { bs, Spectrum::Zero() };
        }

        bs.wo          = Warp::square_to_cosine_hemisphere(sample);
        bs.pdf         = Warp::square_to_cosine_hemisphere_pdf(bs.wo);
        bs.eta         = 1.f;
        bs.sampledType = (uint32_t)EBSDFFlag::DiffuseReflection;

        if (bs.pdf <= 0.)
        {
            return { bs, Spectrum::Zero() };
        }

        auto value = reflectance->evaluate(si);
        return { bs, value };
    }

    Spectrum DiffuseBSDF::evaluate(const SurfaceInteraction& si, const Eigen::Vector3d& wo) const
    {
        double cos_theta_i = Frame::cos_theta(si.wi),
               cos_theta_o = Frame::cos_theta(wo);

        bool active = cos_theta_i > 0. && cos_theta_o > 0.;
        if (!active)
            return Spectrum::Zero();

        auto value = reflectance->evaluate(si) * 1. / EIGEN_PI * cos_theta_o;
        return value;
    }

    double DiffuseBSDF::pdf(const SurfaceInteraction& si, const Eigen::Vector3d& wo) const
    {
        double cos_theta_i = Frame::cos_theta(si.wi),
               cos_theta_o = Frame::cos_theta(wo);

        double pdf = Warp::square_to_cosine_hemisphere_pdf(wo);

        return (cos_theta_i > 0. && cos_theta_o > 0.) ? pdf : 0.;
    }

    std::pair<Spectrum, double> DiffuseBSDF::evaluate_pdf(const SurfaceInteraction& si,
                                                          const Eigen::Vector3d& wo) const
    {
        double cos_theta_i = Frame::cos_theta(si.wi),
               cos_theta_o = Frame::cos_theta(wo);

        bool active = cos_theta_i > 0. && cos_theta_o > 0.;
        if (!active)
            return { Spectrum::Zero(), 0 };

        auto value = reflectance->evaluate(si) * 1. / EIGEN_PI * cos_theta_o;
        double pdf = Warp::square_to_cosine_hemisphere_pdf(wo);
        return { value, pdf };
    }

    Spectrum DiffuseBSDF::diffuseReflectance(const SurfaceInteraction& si) const
    {
        return reflectance->evaluate(si);
    }
}
