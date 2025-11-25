#pragma once

#include "direction_sample.hpp"
#include "position_sample.hpp"
#include "preliminary_interaction.hpp"
#include "ray_fwd.hpp"
#include "transform.hpp"

#include <vislab/core/data.hpp>
#include <memory>

namespace vislab
{
    class BSDF;
    class Light;

    /**
     * @brief Base class for all shapes.
     */
    class Shape : public Data
    {
        VISLAB_INTERFACE(Shape, Object)

    public:
        /**
         * @brief Constructor.
         */
        Shape();

        /**
         * @brief Sample a point on the surface of this shape.
         * @param sample Uniformly distributed 2D point in [0,1]^2.
         * @return Sampled position.
         */
        virtual PositionSample samplePosition(const Eigen::Vector2d& sample) const = 0;

        /**
         * @brief Query the probability density of SamplePosition() for a particular point on the surface.
         * @param ps Position to sample to evaluate.
         * @return Probability density.
         */
        virtual double pdfPosition(const PositionSample& ps) const = 0;

        /**
         * @brief Sample a direction towards this shape with respect to solid angles measured at a reference position within the scene.
         * @param it Interaction from which the sample is taken.
         * @param sample Uniformly distributed 2D point in [0,1]^2.
         * @return Sampled direction.
         */
        virtual DirectionSample sampleDirection(const Interaction& it, const Eigen::Vector2d& sample) const;

        /**
         * @brief Query the probability density of SampleDirection().
         * @param it Interaction from which the sample is taken.
         * @param ds Sampled direction.
         * @return Probability density.
         */
        virtual double pdfDirection(const Interaction& it, const DirectionSample& ds) const;

        /**
         * @brief Fast ray intersection test. To obtain further information about the intersection, use ComputeSurfaceInteraction().
         * @param ray Ray to test against.
         * @return Intersection result.
         */
        virtual PreliminaryIntersection preliminaryHit(const Ray3d& ray) const = 0;

        /**
         * @brief Fast ray shadow test.
         * @param ray Ray to test against.
         * @return True if there was a hit.
         */
        virtual bool anyHit(const Ray3d& ray) const;

        /**
         * @brief Counts the number of hits.
         * @param ray Ray to test against.
         * @return Number of intersections with the shape.
         */
        int32_t countHits(const Ray3d& ray) const;

        /**
         * @brief Compute and return detailed information related to a surface interaction.
         * @param ray Ray to test against.
         * @param pi Preliminary intersection result.
         * @param flags Information about what to query.
         * @return Surface interaction result.
         */
        virtual SurfaceInteraction computeSurfaceInteraction(const Ray3d& ray, const PreliminaryIntersection& pi, EHitComputeFlag flags = EHitComputeFlag::All) const = 0;

        /**
         * @brief Test for an intersection and return detailed information.
         * @param ray Ray to test against.
         * @param flags Information about what to query.
         * @return Surface interaction result.
         */
        SurfaceInteraction closestHit(const Ray3d& ray, EHitComputeFlag flags = EHitComputeFlag::All) const;

        /**
         * @brief Return the shape's surface area.
         * @return Surface area of the shape.
         */
        virtual double surfaceArea() const = 0;

        /**
         * @brief Bounding box of the object in object-space.
         * @return Object-space bounding box.
         */
        virtual Eigen::AlignedBox3d objectBounds() const = 0;

        /**
         * @brief Bounding box the object in world space.
         * @return World-space bounding box.
         */
        virtual Eigen::AlignedBox3d worldBounds() const;

        /**
         * @brief Rebuilds the ray traversal acceleration data structure.
         */
        virtual void buildAccelerationTree();

        /**
         * @brief Transformation of a shape from object space to world space.
         */
        Transform transform;

        /**
         * @brief Bidirectional scattering distribution function.
         */
        std::shared_ptr<BSDF> bsdf;

        /**
         * @brief Bidirectional scattering distribution function.
         */
        std::shared_ptr<Light> light;

        /**
         * @brief Reverses the orientation of the normals.
         */
        bool reverseOrientation;
    };
}
