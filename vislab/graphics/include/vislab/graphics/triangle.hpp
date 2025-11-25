#pragma once

#include "shape.hpp"

#include <vislab/core/array_fwd.hpp>

namespace vislab
{
    class Mesh;

    /**
     * @brief Class that holds a triangle mesh.
     */
    class Triangle : public Shape
    {
        VISLAB_OBJECT(Triangle, Shape)

    public:
        /**
         * @brief Constructor.
         */
        Triangle();

        /**
         * @brief Sample a point on the surface of this shape.
         * @param sample Uniformly distributed 2D point in [0,1]^2.
         * @return Sampled position.
         */
        PositionSample samplePosition(const Eigen::Vector2d& sample) const override;

        /**
         * @brief Query the probability density of SamplePosition() for a particular point on the surface.
         * @param ps Position to sample to evaluate.
         * @return Probability density.
         */
        double pdfPosition(const PositionSample& ps) const override;

        /**
         * @brief Sample a direction towards this shape with respect to solid angles measured at a reference position within the scene.
         * @param it Interaction from which the sample is taken.
         * @param sample Uniformly distributed 2D point in [0,1]^2.
         * @return Sampled direction.
         */
        DirectionSample sampleDirection(const Interaction& it, const Eigen::Vector2d& sample) const override;

        /**
         * @brief Query the probability density of SampleDirection().
         * @param it Interaction from which the sample is taken.
         * @param ds Sampled direction.
         * @return Probability density.
         */
        double pdfDirection(const Interaction& it, const DirectionSample& ds) const override;

        /**
         * @brief Fast ray intersection test. To obtain further information about the intersection, use ComputeSurfaceInteraction().
         * @param ray Ray to test against.
         * @return Intersection result.
         */
        PreliminaryIntersection preliminaryHit(const Ray3d& ray) const override;

        /**
         * @brief Compute and return detailed information related to a surface interaction.
         * @param ray Ray to test against.
         * @param pi Preliminary intersection result.
         * @param flags Information about what to query.
         * @return Surface interaction result.
         */
        SurfaceInteraction computeSurfaceInteraction(const Ray3d& ray, const PreliminaryIntersection& pi, EHitComputeFlag flags = EHitComputeFlag::All) const override;

        /**
         * @brief Return the shape's surface area.
         * @return Surface area of the shape.
         */
        double surfaceArea() const override;

        /**
         * @brief Bounding box of the object in object-space.
         * @return Object-space bounding box.
         */
        Eigen::AlignedBox3d objectBounds() const override;

        /**
         * @brief Rebuilds the ray traversal acceleration data structure.
         */
        void buildAccelerationTree() override;

        /**
         * @brief Mesh
         */
        std::shared_ptr<Mesh> mesh;

    private:
        /**
         * @brief Updates internal variables.
         */
        void update();

        /**
         * @brief Inverse surface area of the sphere.
         */
        double mInverseSurfaceArea;
    };
}
