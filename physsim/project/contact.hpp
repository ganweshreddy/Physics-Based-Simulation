#pragma once

#include <Eigen/Eigen>

namespace physsim
{
    class RigidBody;

    /**
     * @brief Enumeration of contact types.
     */
    enum class ContactType
    {
        Vertex_Face,
        Edge_Edge,
        None
    };

    /**
     * @brief Holds information about a contact.
     */
    struct Contact
    {
        Contact()
            : a(NULL)
            , b(NULL)
            , n(0, 0, 0)
            , p(0, 0, 0)
            , ea(0, 0, 0)
            , eb(0, 0, 0)
            , depth(0)
            , type(ContactType::Vertex_Face)
        {
        }

        /**
         * @brief First object (vertex).
         */
        RigidBody* a;

        /**
         * @brief Second object (face).
         */
        RigidBody* b;

        /**
         * @brief Normal on the face.
         */
        Eigen::Vector3d n;

        /**
         * @brief World-space contact location.
         */
        Eigen::Vector3d p;

        /**
         * @brief Edge direction for A for edge-edge collisions.
         */
        Eigen::Vector3d ea;

        /**
         * @brief Edge direction for B for edge-edge collisions.
         */
        Eigen::Vector3d eb;

        /**
         * @brief Penetration depth.
         */
        double depth;

        /**
         * @brief Type of contact.
         */
        ContactType type;
    };
}
