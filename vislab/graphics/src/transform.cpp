#include <vislab/graphics/transform.hpp>

namespace vislab
{
    Transform::Transform()
        : mMatrix(Eigen::Matrix4d::Identity())
        , mMatrixInverse(Eigen::Matrix4d::Identity())
    {
    }
    Transform::Transform(const Eigen::Matrix4d& matrix)
        : mMatrix(matrix)
        , mMatrixInverse(matrix.inverse())
    {
    }

    const Eigen::Matrix4d& Transform::getMatrix() const { return mMatrix; }
    const Eigen::Matrix4d& Transform::getMatrixInverse() const { return mMatrixInverse; }

    void Transform::setMatrix(const Eigen::Matrix4d& matrix)
    {
        mMatrix        = matrix;
        mMatrixInverse = matrix.inverse();
        matrixChanged.notify(this, &mMatrix);
    }

    Eigen::Vector3d Transform::transformPoint(const Eigen::Vector3d& point) const
    {
        Eigen::Vector4d p4 = mMatrix * Eigen::Vector4d(point.x(), point.y(), point.z(), 1.);
        return p4.xyz() / p4.w();
    }

    Eigen::Vector3d Transform::transformVector(const Eigen::Vector3d& direction) const
    {
        return mMatrix.block(0, 0, 3, 3) * direction;
    }

    Eigen::Vector3d Transform::transformNormal(const Eigen::Vector3d& normal) const
    {
        return (mMatrixInverse.block(0, 0, 3, 3).transpose() * normal).stableNormalized();
    }

    Ray3d Transform::transformRay(const Ray3d& ray) const
    {
        Ray3d result     = ray;
        result.origin    = transformPoint(ray.origin);
        result.direction = transformVector(ray.direction);
        return result;
    }

    Eigen::Vector3d Transform::transformPointInverse(const Eigen::Vector3d& point) const
    {
        Eigen::Vector4d p4 = mMatrixInverse * Eigen::Vector4d(point.x(), point.y(), point.z(), 1.);
        return p4.xyz() / p4.w();
    }

    Eigen::Vector3d Transform::transformVectorInverse(const Eigen::Vector3d& direction) const
    {
        return mMatrixInverse.block(0, 0, 3, 3) * direction;
    }

    Ray3d Transform::transformRayInverse(const Ray3d& ray) const
    {
        Ray3d result     = ray;
        result.origin    = transformPointInverse(ray.origin);
        result.direction = transformVectorInverse(ray.direction);
        return result;
    }
}
