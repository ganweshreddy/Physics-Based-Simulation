#pragma once

#include "bsdf_sample.hpp"
#include "spectrum.hpp"
#include "surface_interaction.hpp"

#include <vislab/core/data.hpp>

#include <Eigen/Eigen>

namespace vislab
{
    /**
     * @brief Base class for bidirectional scattering distribution functions.
     */
    class BSDF : public Data
    {
        VISLAB_INTERFACE(BSDF, Data)

    public:
        /**
         * @brief Constructor.
         * @param flags Flags that characterize the requirements of the BSDF.
         */
        BSDF(uint32_t flags);

        /**
         * @brief Importance samples the BSDF model.
         * @param si Surface interaction that describes where to sample the BSDF.
         * @param sample Uniformly distributed sample in [0,1]^2.
         * @return Sampling record and BSDF value (divided by PDF for non-delta BSDFs).
         */
        virtual std::pair<BSDFSample, Spectrum> sample(const SurfaceInteraction& si, const Eigen::Vector2d& sample) const = 0;

        /**
         * @brief Evaluates the BSDF for a given outgoing direction.
         * @param si Surface interaction that describes where to evaluate the BSDF.
         * @param wo Outgoing direction in local coordinates.
         * @return BSDF evaluation.
         */
        virtual Spectrum evaluate(const SurfaceInteraction& si, const Eigen::Vector3d& wo) const = 0;

        /**
         * @brief Probability per unit solid angle for sampling a given direction.
         * @param si Surface interaction that describes where to evaluate the BSDF.
         * @param wo Outgoing direction in local coordinates.
         * @return Probability per unit solid angle.
         */
        virtual double pdf(const SurfaceInteraction& si, const Eigen::Vector3d& wo) const = 0;

        /**
         * @brief Jointly evaluates the BSDF and the probability density.
         * @param si Surface interaction that describes where to evaluate the BSDF.
         * @param wo Outgoing direction in local coordinates.
         * @return Evaluated BSDF and probability density.
         */
        virtual std::pair<Spectrum, double> evaluate_pdf(const SurfaceInteraction& si, const Eigen::Vector3d& wo) const;

        /**
         * @brief Jointly evaluates the BSDF, the probability density, and importance samples the BSDF.
         * @param si Surface interaction that describes where to evaluate the BSDF.
         * @param wo Outgoing direction in local coordinates.
         * @param sample Uniformly distributed sample in [0,1]^2.
         * @return Evaluated BSDF, probability density, importance sampling result.
         */
        virtual std::tuple<Spectrum, double, BSDFSample, Spectrum> evaluate_pdf_sample(const SurfaceInteraction& si, const Eigen::Vector3d& wo, const Eigen::Vector2d& sample) const;

        /**
         * @brief Evaluates the diffuse reflectance.
         * @param si Surface interaction that describes where to evaluate the BSDF.
         * @return Diffuse reflectance.
         */
        virtual Spectrum diffuseReflectance(const SurfaceInteraction& si) const;

        /**
         * @brief Flags that characterize the requirements for evaluating this BSDF.
         * @return Flags.
         */
        uint32_t flags() const;

        /**
         * @brief Does the implementation require access to texture-space differentials?
         * @return True if differentials are needed for the BSDF evaluation.
         */
        bool needsDifferentials() const;

        /**
         * @brief Tests if a flag is set.
         * @param flag Flag to test.
         * @return True, if the flag is set.
         */
        bool hasFlag(const EBSDFFlag& flag) const;

    protected:
        /**
         * @brief Flags that characterize the requirements for evaluating this BSDF.
         */
        const uint32_t mFlags;
    };
}
