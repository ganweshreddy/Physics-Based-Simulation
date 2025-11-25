#pragma once

#include "light.hpp"

namespace vislab
{
    /**
     * @brief Point light source.
     */
    class PointLight : public Light
    {
        VISLAB_OBJECT(PointLight, Light)

    public:
        /**
         * @brief Constructor.
         */
        PointLight();

        /**
         * @brief Generates a ray in world space.
         * @param sample Sample location in [0,1]^2.
         * @param sample2 Another sample location in [0,1]^2.
         * @return Ray that was sampled and spectrum.
         */
        std::pair<Ray3d, Spectrum> sampleRay(const Eigen::Vector2d& sample, const Eigen::Vector2d& sample2) const override;

        /**
         * @brief Given a reference point in the scene, sample a direction from the reference point towards the endpoint (ideally proportional to the emission/sensitivity profile).
         * @param it Interaction to sample the direction from.
         * @param sample Sample location in [0,1]^2.
         * @return Information about the sampled direction.
         */
        std::pair<DirectionSample, Spectrum> sampleDirection(const Interaction& it, const Eigen::Vector2d& sample) const override;

        /**
         * @brief Evaluate the probability density of a direct sampling method implemented by the sampleDirection() method.
         * @param it Interaction to sample the direction from.
         * @param ds Information about the sampled direction.
         * @return Probability density.
         */
        double pdfDirection(const Interaction& it, const DirectionSample& ds) const override;

        /**
         * @brief Given a ray-surface intersection, return the emitted radiance or importance traveling along the reverse direction.
         * @param si Surface interaction that defines where to evaluate.
         * @return Evaluated radiance or importance.
         */
        Spectrum evaluate(const SurfaceInteraction& si) const override;

        /**
         * @brief Flags that characterize the requirements for evaluating this light source.
         * @return Flags.
         */
        uint32_t flags() const override;

        /**
         * @brief Radiant intensity "I", i.e., radiant flux per unit solid angle.
         */
        Eigen::Vector3d intensity;
    };
}
