#include "rigid_body.hpp"

#include <vislab/graphics/shape.hpp>
#include <vislab/graphics/transform.hpp>
#include <vislab/graphics/triangle.hpp>
#include <vislab/graphics/mesh.hpp>
#include <vislab/core/array.hpp>

namespace physsim
{
    RigidBody::RigidBody()
        : mPosition(Eigen::Vector3d::Zero())
        , mRotation(Eigen::Quaterniond::Identity())
        , mLinearVelocity(Eigen::Vector3d::Zero())
        , mAngularVelocity(Eigen::Vector3d::Zero())
        , mInertiaBody(Eigen::Matrix3d::Identity())
        , mInertiaBodyInverse(Eigen::Matrix3d::Identity())
        , mMass(1)
        , mMassInverse(1)
        , mForce(Eigen::Vector3d::Zero())
        , mTorque(Eigen::Vector3d::Zero())
        , mType(EType::Dynamic)
        , mShape(nullptr)
        , mScale(1)
    {
    }

    const Eigen::Vector3d& RigidBody::position() const { return mPosition; }

    void RigidBody::setPosition(const Eigen::Vector3d& x)
    {
        mPosition = x;
        updateTransform();
    }

    const Eigen::Quaterniond& RigidBody::rotation() const { return mRotation; }

    void RigidBody::setRotation(const Eigen::Quaterniond& q)
    {
        mRotation = q;
        updateTransform();
    }

    Eigen::Vector3d RigidBody::linearMomentum() const
    {
        return mMass * mLinearVelocity;
    }

    Eigen::Vector3d RigidBody::angularMomentum() const
    {
        return inertiaWorld() * mAngularVelocity;
    }

    const Eigen::Matrix3d& RigidBody::inertiaBody() const { return mInertiaBody; }

    void RigidBody::setInertiaBody(const Eigen::Matrix3d& Ib)
    {
        mInertiaBody        = Ib;
        mInertiaBodyInverse = Ib.inverse();
    }

    const Eigen::Matrix3d& RigidBody::inertiaBodyInverse() const { return mInertiaBodyInverse; }

    const double& RigidBody::mass() const { return mMass; }

    void RigidBody::setMass(const double& m)
    {
        mMass        = m;
        mMassInverse = 1. / m;
    }

    const double& RigidBody::massInverse() const { return mMassInverse; }

    const Eigen::Vector3d& RigidBody::force() const { return mForce; }

    const Eigen::Vector3d& RigidBody::torque() const { return mTorque; }

    const RigidBody::EType& RigidBody::type() const { return mType; }

    void RigidBody::setType(const EType& type) { mType = type; }

    Eigen::Matrix3d RigidBody::inertiaWorld() const
    {
        return mRotation * mInertiaBody * mRotation.inverse();
    }

    Eigen::Matrix3d RigidBody::inertiaWorldInverse() const
    {
        return mRotation * mInertiaBodyInverse * mRotation.inverse();
    }

    const Eigen::Vector3d& RigidBody::linearVelocity() const
    {
        return mLinearVelocity;
    }

    void RigidBody::setLinearVelocity(const Eigen::Vector3d& v)
    {
        if (mType != EType::Dynamic)
            return;
        moldLinearVelocity = std::exchange(mLinearVelocity, v);
    }

    double RigidBody::veldiff()
    {
        return (mLinearVelocity - moldLinearVelocity).norm();
    }

    const Eigen::Vector3d& RigidBody::angularVelocity() const
    {
        return mAngularVelocity;
    }

    void RigidBody::setAngularVelocity(const Eigen::Vector3d& w)
    {
        if (mType != EType::Dynamic)
            return;
        mAngularVelocity = w;
    }

    Eigen::Vector3d RigidBody::velocity(const Eigen::Vector3d& point) const
    {
        return linearVelocity() +
               angularVelocity().cross(point - position());
    }

    void RigidBody::applyForceToCenterOfMass(const Eigen::Vector3d& f)
    {
        if (mType != EType::Dynamic)
            return;
        mForce += f;
    }

    void RigidBody::applyForce(const Eigen::Vector3d& f, const Eigen::Vector3d& p)
    {
        if (mType != EType::Dynamic)
            return;
        mForce += f;
        mTorque += (p - mPosition).cross(f);
    }

    void RigidBody::applyTorque(const Eigen::Vector3d& t)
    {
        if (mType != EType::Dynamic)
            return;
        mTorque += t;
    }

    void RigidBody::resetForce()
    {
        if (mType != EType::Dynamic)
            return;
        mForce = Eigen::Vector3d::Zero();
    }

    void RigidBody::resetTorque()
    {
        if (mType != EType::Dynamic)
            return;
        mTorque = Eigen::Vector3d::Zero();
    }

    std::shared_ptr<const vislab::Shape> RigidBody::shape() const
    {
        return mShape;
    }

    void RigidBody::setShape(std::shared_ptr<vislab::Shape> shape)
    {
        mShape = shape;
        updateTransform();
    }

    const double& RigidBody::scale() const
    {
        return mScale;
    }

    void RigidBody::setScale(const double& scale)
    {
        mScale = scale;
        updateTransform();
    }

    void RigidBody::updateTransform()
    {
        if (mShape)
        {
            Eigen::Matrix4d transform   = Eigen::Matrix4d::Identity();
            transform.block(0, 0, 3, 3) = mRotation.toRotationMatrix();
            transform.block(0, 3, 3, 1) = mPosition;
            mShape->transform.setMatrix(transform * Eigen::Matrix4d::scale({ mScale, mScale, mScale }));
        }
    }
}
