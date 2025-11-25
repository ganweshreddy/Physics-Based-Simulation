#pragma once

#include "shape.hpp"

namespace vislab
{
    /**
     * @brief Class that defines a rectangle spanning (-1,-1,0) to (1,1,0). To modify its position update the transformation.
     */
    class Rectangle : public Shape
    {
        VISLAB_OBJECT(Rectangle, Shape)

    public:
        /**
         * @brief Constructor.
         */
        Rectangle();

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

    private:
        /**
         * @brief Updates internal variables.
         */
        void update();

        /**
         * @brief Inverse surface area of the rectangle.
         */
        double mInverseSurfaceArea;

        /**
         * @brief Local coordinate frame.
         */
        Frame mFrame;
    };
}
