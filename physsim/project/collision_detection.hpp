#pragma once

#include "contact.hpp"
#include "rigid_body.hpp"

#include <vislab/core/array.hpp>
#include <vislab/graphics/mesh.hpp>
#include <vislab/graphics/triangle.hpp>

#include <Eigen/Eigen>
#include <set>
#include <utility>
#include <vector>

namespace physsim
{
    struct TransformedMesh;

    /**
     * @brief Enumeration of broad phase collision detection methods.
     */
    enum class EBroadPhaseMethod
    {
        None,
        AABB,
        SweepAndPrune
    };

    /**
     * @brief Enumeration of narrow phase collision detection methods.
     */
    enum class ENarrowPhaseMethod
    {
        Exhaustive,
        GilbertJohnsonKeerthi
    };

    /**
     * @brief Class to detect and resolve collisions between rigid bodies.
     */
    class CollisionDetection
    {
    public:
        CollisionDetection(const std::vector<std::shared_ptr<RigidBody>>& objects);

        /**
         * @brief Performs a collision detection.
         * @param broadPhaseMethod Broad phase detection method.
         * @param narrowPhaseMethod Narrow phase detection method.
         * @param eps Impulse epsilon, typically in [0,1]. A value of 1 is energy conserving. Everything below adds dampening.
         * @param stepSize Numerical integration step size.
         */
        void computeCollisionDetection(EBroadPhaseMethod broadPhaseMethod, ENarrowPhaseMethod narrowPhaseMethod, double eps, double stepSize);

    private:
        /**
         * @brief Perform the broad phase collision detection.
         * @param broadPhaseMethod Broad phase method to apply.
         */
        void computeBroadPhase(EBroadPhaseMethod broadPhaseMethod);

        /**
         * @brief Perform narrow phase collision detection.
         * @param narrowPhaseMethod Narrow phase method to apply.
         */
        void computeNarrowPhase(ENarrowPhaseMethod narrowPhaseMethod);

        // test if ray from start towards end intersects object with vertices V and faces F
        ContactType isColliding(const vislab::Ray3d& ray, const vislab::Shape* shape);

        /**
         * @brief Loop over edges of faces and compute closest to given edge.
         * @param start Start of given edge.
         * @param end End of given edge.
         * @param V Set of vertices.
         * @param F Set of face indices.
         * @return Contact information.
         */
        Contact findEdgeEdgeCollision(const Eigen::Vector3d& start, const Eigen::Vector3d& end, const TransformedMesh& mesh);

        /**
         * @brief Apply impulses for the found contacts.
         * @param eps Impulse epsilon.
         * @param stepSize Numerical integration step size.
         */
        void applyImpulse(double eps, double stepSize);

        /**
         * @brief Clear all data structures.
         */
        void clearDataStructures();

        /**
         * @brief All rigid body objects in the scene.
         */
        std::vector<std::shared_ptr<RigidBody>> mObjects;

        /**
         * @brief Result of broadphase, pairs of objects with possible collisions.
         */
        std::vector<std::pair<size_t, size_t>> mOverlappingBodys;

        /**
         * @brief Set of vertex indices that penetrate a face, used to avoid duplicates.
         */
        std::set<int> mPenetratingVertices;

        /**
         * @brief Set of pairs of vertex indices that represent a penetrating edge, used to avoid duplicates.
         */
        std::set<std::pair<int, int>> mPenetratingEdges;

        /**
         * @brief Computed contact points.
         */
        std::vector<Contact> mContacts;
    };
}
