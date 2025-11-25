#pragma once

#include "bsdf.hpp"

#include <Eigen/Eigen>
#include <memory>

namespace vislab
{
    class Texture;

    /**
     * @brief A BSDF for diffuse Lambertian shading.
     */
    class DiffuseBSDF : public BSDF
    {
        VISLAB_OBJECT(DiffuseBSDF, BSDF)

    public:
        /**
         * @brief Constructor.
         */
        DiffuseBSDF();

        /**
         * @brief Constructor with initial reflectance.
         */
        DiffuseBSDF(std::shared_ptr<Texture> reflectance);

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
         * @brief Jointly evaluates the BSDF and the probability density.
         * @param si Surface interaction describes the sampling location.
         * @param wo Outgoing direction.
         * @return BSDF and probability density.
         */
        std::pair<Spectrum, double> evaluate_pdf(const SurfaceInteraction& si, const Eigen::Vector3d& wo) const override;

        /**
         * @brief Evaluates the diffuse reflectance.
         * @param si Surface interaction describes the sampling location.
         * @return Approximation as diffuse material.
         */
        Spectrum diffuseReflectance(const SurfaceInteraction& si) const override;

        /**
         * @brief Reflectance of the material.
         */
        std::shared_ptr<Texture> reflectance;
    };
}
