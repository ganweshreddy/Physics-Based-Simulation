#pragma once

#include "contact.hpp"

namespace physsim
{
    struct TransformedMesh;

    /**
     * @brief Implementation of the GJK algorithm.
     */
    class GilbertJohnsonKeerthi
    {
    public:
        /**
         * @brief Determines whether there is a collision between two sets of vertices.
         * @param A Vertex set A.
         * @param B Vertex set B.
         * @return Contact information.
         */
        static Contact findCollisionGJK(const TransformedMesh& A, const TransformedMesh& B);

    private:
        /**
         * @brief Stores the support function results of both shapes and the Minkowski differences.
         */
        struct SupportVec
        {
            /**
             * @brief Empty constructor.
             */
            SupportVec()
                : A(0, 0, 0)
                , B(0, 0, 0)
                , P(0, 0, 0)
            {
            }

            /**
             * @brief Constructor from two points and the difference of the points.
             * @param a Furthest point on A in direction dir.
             * @param b Furthest point on B in direction -dir.
             * @param p Minkowski difference.
             */
            SupportVec(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& p)
                : A(a)
                , B(b)
                , P(p)
            {
            }

            /**
             * @brief Support point on A in dir.
             */
            Eigen::Vector3d A;

            /**
             * @brief Support point on B in -dir.
             */
            Eigen::Vector3d B;

            /**
             * @brief Minkowski difference of A and B.
             */
            Eigen::Vector3d P;
        };

        /**
         * @brief Locates the vertex from a set of vertices V furthest into a direction d.
         * @param V Matrix that contains the vertices.
         * @param d Direction vector.
         * @return Coordinate of the farthest vertex in direction.
         */
        static Eigen::Vector3d support(const TransformedMesh& mesh, const Eigen::Vector3d& d);

        /**
         * @brief Locates the vertex on the Minkowski difference of sets of vertices A and B furthest into a direction d
         * @param A First vertex set.
         * @param B Second vertex set.
         * @param d Direction to compute difference for.
         * @return Support vector.
         */
        static SupportVec support(const TransformedMesh& A, const TransformedMesh& B, const Eigen::Vector3d& d);

        /**
         * @brief Calculates the closest simplex of a simplex S=(a,b,c,d), where dim is the number of dimensions and a is the last added vertex. Based on: https://github.com/kevinmoran/GJK
         * @param a First point of simplex.
         * @param b Second point of simplex.
         * @param c Third point of simplex.
         * @param d Fourth point of simplex.
         * @param dim Number of dimensions.
         * @param search_dir Search direction.
         * @return True if origin is contained in simplex.
         */
        static bool nearestSimplex(SupportVec& a, SupportVec& b, SupportVec& c, SupportVec& d, int& dim, Eigen::Vector3d& search_dir);

        /**
         * @brief Expanding Polytope Algorithm. Based on: https://github.com/kevinmoran/GJK
         * @param a First point of simplex.
         * @param b Second point of simplex.
         * @param c Third point of simplex.
         * @param d Fourth point of simplex.
         * @param A Vertex set A.
         * @param B Vertex set B.
         * @return Contact information.
         */
        static Contact EPA(const SupportVec& a, const SupportVec& b, const SupportVec& c, const SupportVec& d, const TransformedMesh& A, const TransformedMesh& B);
    };
}
