#pragma once

#include "bsdf.hpp"

#include <Eigen/Eigen>
#include <memory>

namespace vislab
{
    class Texture;

    /**
     * @brief A BSDF for a smooth dielectric material.
     */
    class DielectricBSDF : public BSDF
    {
        VISLAB_OBJECT(DielectricBSDF, BSDF)

    public:
        /**
         * @brief Constructor.
         */
        DielectricBSDF();

        /**
         * @brief Importance sample the BSDF model.
         * @param si Surface interaction describes the sampling location.
         * @param sample Uniformly distributed random number in [0,1]^2.
         * @return Sampling record and BSDF divided by the probability.
         */
        std::pair<BSDFSample, Spectrum> sample(const SurfaceInteraction& si, const Eigen::Vector2d& sample) const override;

        /**
         * @brief Evaluation of the BSDF for a given outgoing direction.
         * @param si Surface interaction describes the sampling location.
         * @param wo Outgoing direction.
         * @return Evaluated BSDF.
         */
        Spectrum evaluate(const SurfaceInteraction& si, const Eigen::Vector3d& wo) const override;

        /**
         * @brief Probability per unit solid angle for sampling a given direction.
         * @param si Surface interaction describes the sampling location.
         * @param wo Outgoing direction.
         * @return Probability density.
         */
        double pdf(const SurfaceInteraction& si, const Eigen::Vector3d& wo) const override;

        /**
         * @brief Interior index of reflection.
         */
        double interior_ior;

        /**
         * @brief Exterios index of reflection.
         */
        double exterior_ior;

        /**
         * @brief Relative index of refraction from the exterior to the interior.
         */
        double eta;

        /**
         * @brief Specular reflectance of the material.
         */
        std::shared_ptr<Texture> specularReflectance;

        /**
         * @brief Specular transmittance of the material.
         */
        std::shared_ptr<Texture> specularTransmittance;
    };
}
