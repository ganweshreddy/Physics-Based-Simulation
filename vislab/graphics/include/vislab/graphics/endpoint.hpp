#pragma once

#include "direction_sample.hpp"
#include "ray_fwd.hpp"
#include "spectrum.hpp"

#include <vislab/core/data.hpp>

#include <Eigen/Eigen>
#include <memory>

namespace vislab
{
    class Medium;
    class Shape;
    class Interaction;

    /**
     * @brief Base class for emitters and sensors.
     */
    class Endpoint : public Data
    {
        VISLAB_INTERFACE(Endpoint, Data)

    public:
        /**
         * @brief Generates a ray in world space.
         * @param sample Sample location in [0,1]^2.
         * @param sample2 Another sample location in [0,1]^2.
         * @return Ray that was sampled and spectrum.
         */
        virtual std::pair<Ray3d, Spectrum> sampleRay(const Eigen::Vector2d& sample, const Eigen::Vector2d& sample2) const = 0;

        /**
         * @brief Given a reference point in the scene, sample a direction from the reference point towards the endpoint (ideally proportional to the emission/sensitivity profile).
         * @param it Interaction to sample the direction from.
         * @param sample Sample location in [0,1]^2.
         * @return Information about the sampled direction.
         */
        virtual std::pair<DirectionSample, Spectrum> sampleDirection(const Interaction& it, const Eigen::Vector2d& sample) const;

        /**
         * @brief Evaluate the probability density of a direct sampling method implemented by the sampleDirection() method.
         * @param it Interaction to sample the direction from.
         * @param ds Information about the sampled direction.
         * @return Probability density.
         */
        virtual double pdfDirection(const Interaction& it, const DirectionSample& ds) const;

        /**
         * @brief Given a ray-surface intersection, return the emitted radiance or importance traveling along the reverse direction.
         * @param si Surface interaction that defines where to evaluate.
         * @return Evaluated radiance or importance.
         */
        virtual Spectrum evaluate(const SurfaceInteraction& si) const;

        /**
         * @brief Shape, to which the end point is currently attached.
         */
        std::shared_ptr<Shape> shape;

        /**
         * @brief Medium that surrounds the end point.
         */
        std::shared_ptr<Medium> medium;
    };
}
