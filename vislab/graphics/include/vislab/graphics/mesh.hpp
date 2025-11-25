#pragma once

#include "bounding_volume_hierarchy.hpp"
#include "shape.hpp"

#include <vislab/core/array_fwd.hpp>

#include <vislab/core/event.hpp>

namespace vislab
{
    /**
     * @brief Class that holds a mesh.
     */
    class Mesh : public Data
    {
        VISLAB_OBJECT(Mesh, Data)

    public:
        /**
         * @brief Constructor.
         */
        Mesh();

        /**
         * @brief Sample a point on the surface of this shape.
         * @param sample Uniformly distributed 2D point in [0,1]^2.
         * @return Sampled position.
         */
        PositionSample samplePosition(const Eigen::Vector2d& sample) const;

        /**
         * @brief Query the probability density of SamplePosition() for a particular point on the surface.
         * @param ps Position to sample to evaluate.
         * @return Probability density.
         */
        double pdfPosition(const PositionSample& ps) const;

        /**
         * @brief Sample a direction towards this shape with respect to solid angles measured at a reference position within the scene.
         * @param it Interaction from which the sample is taken.
         * @param sample Uniformly distributed 2D point in [0,1]^2.
         * @return Sampled direction.
         */
        DirectionSample sampleDirection(const Interaction& it, const Eigen::Vector2d& sample) const;

        /**
         * @brief Query the probability density of SampleDirection().
         * @param it Interaction from which the sample is taken.
         * @param ds Sampled direction.
         * @return Probability density.
         */
        double pdfDirection(const Interaction& it, const DirectionSample& ds) const;

        /**
         * @brief Fast ray intersection test. To obtain further information about the intersection, use ComputeSurfaceInteraction().
         * @param ray Ray to test against.
         * @return Intersection result.
         */
        PreliminaryIntersection preliminaryHit(const Ray3d& ray) const;

        /**
         * @brief Compute and return detailed information related to a surface interaction.
         * @param ray Ray to test against.
         * @param pi Preliminary intersection result.
         * @param flags Information about what to query.
         * @return Surface interaction result.
         */
        SurfaceInteraction computeSurfaceInteraction(const Ray3d& ray, const PreliminaryIntersection& pi, EHitComputeFlag flags = EHitComputeFlag::All) const;

        /**
         * @brief Return the shape's surface area.
         * @return Surface area of the shape.
         */
        double surfaceArea() const;

        /**
         * @brief Bounding box of the object in object-space.
         * @return Object-space bounding box.
         */
        Eigen::AlignedBox3d objectBounds() const;

        /**
         * @brief Performs a ray triangle intersection.
         * @param faceIndex Triangle to perform intersection with.
         * @param ray Ray to intersect with.
         * @return Intersection information.
         */
        PreliminaryIntersection rayTriangleIntersection(const uint32_t& faceIndex, const Ray3d& ray) const;

        /**
         * @brief Recomputes the vertex normals from positions.
         */
        void recomputeVertexNormals();

        /**
         * @brief Mesh vertex positions in object space.
         */
        std::shared_ptr<Array3f> positions;

        /**
         * @brief Mesh vertex normals.
         */
        std::shared_ptr<Array3f> normals;

        /**
         * @brief Mesh texture coordinates.
         */
        std::shared_ptr<Array2f> texCoords;

        /**
         * @brief Index buffer.
         */
        std::shared_ptr<Array3ui> indices;

        /**
         * @brief Event that can be raised when the positions changed.
         */
        TEvent<Mesh, void> positionsChanged;

        /**
         * @brief Rebuilds the ray traversal acceleration data structure.
         */
        void buildAccelerationTree();

    private:
        /**
         * @brief Triangle leaf element in the bounding volume hierarchy.
         */
        class Triangle
        {
        public:
            /**
             * @brief Constructor.
             */
            Triangle();

            /**
             * @brief World bounding box of the triangle.
             * @return Bounding box in world space.
             */
            Eigen::AlignedBox3d worldBounds() const;

            /**
             * @brief Performs a ray intersection test against the triangle.
             * @param ray Ray to test against.
             * @return Intersection information.
             */
            PreliminaryIntersection preliminaryHit(const Ray3d& ray) const;

            /**
             * @brief Face index of the triangle in the mesh.
             */
            uint32_t face;

            /**
             * @brief Pointer to the parent mesh.
             */
            Mesh* mesh;

            /**
             * @brief Bounding box of the triangle.
             */
            Eigen::AlignedBox3d worldBoundingBox;
        };

        /**
         * @brief Bounding volume hierarchy for fast ray traversal.
         */
        std::shared_ptr<BoundingVolumeHierarchy<Triangle, 3>> mBoundingVolumeHierarchy;

        /**
         * @brief Storage of all the triangles in the bounding volume hierarchy.
         */
        std::vector<Triangle> mTriangles;
    };
}
