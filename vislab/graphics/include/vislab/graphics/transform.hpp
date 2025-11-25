#pragma once

#include "ray.hpp"

#include <vislab/core/data.hpp>
#include <vislab/core/event.hpp>

#include <Eigen/Eigen>

namespace vislab
{
    /**
     * @brief Stores a transformation.
     */
    class Transform : public Data
    {
        VISLAB_OBJECT(Transform, Data)

    public:
        /**
         * @brief Constructor with an identity transformation.
         */
        Transform();

        /**
         * @brief Constructor with a given initial matrix.
         * @param matrix Initial matrix.
         */
        Transform(const Eigen::Matrix4d& matrix);

        /**
         * @brief Gets the transformation matrix.
         * @return Transformation matrix.
         */
        const Eigen::Matrix4d& getMatrix() const;

        /**
         * @brief Inverse of the transformation matrix.
         * @return Inverse matrix.
         */
        const Eigen::Matrix4d& getMatrixInverse() const;

        /**
         * @brief Sets the world transformation matrix.
         * @param matrix Transformation matrix to set.
         */
        void setMatrix(const Eigen::Matrix4d& matrix);

        /**
         * @brief Event that is raised when the matrix changed.
         */
        TEvent<Transform, Eigen::Matrix4d> matrixChanged;

        /**
         * @brief Affine transformation of a 3D coordinate.
         * @param point Point to transform.
         * @return Transformed point.
         */
        Eigen::Vector3d transformPoint(const Eigen::Vector3d& point) const;

        /**
         * @brief Transforms a direction vector.
         * @param direction Direction vector to transform.
         * @return Transformed direction vector.
         */
        Eigen::Vector3d transformVector(const Eigen::Vector3d& direction) const;

        /**
         * @brief Transforms a normal vector.
         * @param normal Normal vector to transform.
         * @return Transformed normal vector.
         */
        Eigen::Vector3d transformNormal(const Eigen::Vector3d& normal) const;

        /**
         * @brief Transforms a ray by the matrix.
         * @param ray Ray to transform.
         * @return Transformed ray.
         */
        Ray3d transformRay(const Ray3d& ray) const;

        /**
         * @brief Affine inverse transformation of a 3D coordinate.
         * @param point Point to transform.
         * @return Transformed point.
         */
        Eigen::Vector3d transformPointInverse(const Eigen::Vector3d& point) const;

        /**
         * @brief Inverse transform of a direction vector.
         * @param direction Direction vector to transform.
         * @return Transformed direction vector.
         */
        Eigen::Vector3d transformVectorInverse(const Eigen::Vector3d& direction) const;

        /**
         * @brief Inverse transform of a ray by the matrix.
         * @param ray Ray to transform.
         * @return Transformed ray.
         */
        Ray3d transformRayInverse(const Ray3d& ray) const;

    private:
        /**
         * @brief Transformation matrix.
         */
        Eigen::Matrix4d mMatrix;

        /**
         * @brief Inverse of transformation matrix.
         */
        Eigen::Matrix4d mMatrixInverse;
    };
}
