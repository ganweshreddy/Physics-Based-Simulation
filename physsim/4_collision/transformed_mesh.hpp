#pragma once

#include <vislab/core/array.hpp>
#include <vislab/graphics/transform.hpp>

#include <Eigen/Eigen>

namespace physsim
{
    /**
     * @brief Helper class that accesses vertices in world space.
     */
    struct TransformedMesh
    {
        /**
         * @brief Constructor.
         * @param t Object to world transformation.
         * @param v Vertices in object space.
         * @param i Triangle indices.
         */
        TransformedMesh(const vislab::Transform& t,
                        const vislab::Array3f& v,
                        const vislab::Array3ui& i)
            : transform(t)
            , vertices(v)
            , indices(i)
        {
        }

        /**
         * @brief Gets a vertex in world coordinates.
         * @param index Index of vertex to get.
         * @return Vertex in world coordinates.
         */
        Eigen::Vector3d getVertex(const Eigen::Index& index) const
        {
            return transform.transformPoint(vertices.getValueDouble(index));
        }

        /**
         * @brief Gets a primitive by its triangle index.
         * @param index Triangle index.
         * @return Triplet of vertex corners.
         */
        typename vislab::Array3ui::ConstElementRef getPrimitive(const Eigen::Index& index) const
        {
            return indices.getValue(index);
        }

        /**
         * @brief Object to world transformation.
         */
        const vislab::Transform& transform;

        /**
         * @brief Vertices in object space.
         */
        const vislab::Array3f& vertices;

        /**
         * @brief Triangle indices.
         */
        const vislab::Array3ui& indices;
    };
}
