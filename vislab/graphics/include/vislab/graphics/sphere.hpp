#pragma once

#include "shape.hpp"

namespace vislab
{
    /**
     * @brief Class that defines a unit sphere at the origin (0,0,0). To modify its radius or center position update the transformation. Note that non-uniform scaling is not allowed.
     */
    class Sphere : public Shape
    {
        VISLAB_OBJECT(Sphere, Shape)

    public:
        /**
         * @brief Constructor.
         */
        Sphere();

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
         * @brief Fast ray shadow test.
         * @param ray Ray to test against.
         * @return True if there was a hit.
         */
        bool anyHit(const Ray3d& ray) const override;

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
         * @brief Gets the radius of the sphere. To change the radius, update the transformation. Note that non-uniform scalings are not allowed.
         * @return Radius of the sphere.
         */
        double radius() const;

        /**
         * @brief Gets the center of the sphere in world coordinates. To change the center, update the transformation.
         * @return Center location of the sphere.
         */
        Eigen::Vector3d center() const;

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
