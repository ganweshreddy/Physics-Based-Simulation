#include <vislab/graphics/projective_camera.hpp>

#include <vislab/graphics/interaction.hpp>
#include <vislab/graphics/ray.hpp>

namespace vislab
{
    ProjectiveCamera::ProjectiveCamera()
        : mProj(Eigen::Matrix4d::Identity())
        , mNear(0.01)
        , mFar(1.0)
        , mWidth(1.0)
        , mHeight(1.0)
    {
    }

    const Eigen::Matrix4d& ProjectiveCamera::getProj() const { return mProj; }
    const double& ProjectiveCamera::getNear() const { return mNear; }
    const double& ProjectiveCamera::getFar() const { return mFar; }
    const double& ProjectiveCamera::getWidth() const { return mWidth; }
    const double& ProjectiveCamera::getHeight() const { return mHeight; }
    double ProjectiveCamera::getAspectRatio() const { return mWidth / mHeight; }

    void ProjectiveCamera::setNear(const double& near)
    {
        mNear = near;
        updateProjMatrix();
    }

    void ProjectiveCamera::setFar(const double& far)
    {
        mFar = far;
        updateProjMatrix();
    }

    void ProjectiveCamera::setWidth(const double& width)
    {
        mWidth = width;
        updateProjMatrix();
    }

    void ProjectiveCamera::setHeight(const double& height)
    {
        mHeight = height;
        updateProjMatrix();
    }

    Eigen::Vector2d ProjectiveCamera::project(const Eigen::Vector3d& coord) const
    {
        Eigen::Vector4d p   = Eigen::Vector4d(coord.x(), coord.y(), coord.z(), 1);
        Eigen::Vector4d npc = mProj * (mView * p);
        npc /= npc.w();
        return Eigen::Vector2d(
            (npc.x() * 0.5 + 0.5) * mWidth,
            (-npc.y() * 0.5 + 0.5) * mHeight);
    }

    Eigen::Vector3d ProjectiveCamera::unproject(const Eigen::Vector3d& coord) const
    {
        const double viewport_minx = 0;
        const double viewport_miny = 0;
        Eigen::Matrix4d matrix     = mProj * mView;
        Eigen::Vector3d v(
            (((coord.x() - viewport_minx) / mWidth) * 2.0) - 1.0,
            -((((coord.y() - viewport_miny) / mHeight) * 2.0) - 1.0),
            (coord.z() - mNear) / (mFar - mNear));
        Eigen::Vector4d vv = matrix.inverse() * Eigen::Vector4d(v.x(), v.y(), v.z(), 1);
        return vv.xyz() / vv.w();
    }

    std::pair<Ray3d, Spectrum> ProjectiveCamera::sampleRay(const Eigen::Vector2d& sample, const Eigen::Vector2d& sample2) const
    {
        Eigen::Vector3d pixelNear(mWidth * sample.x(), mHeight * (1 - sample.y()), mNear);
        Eigen::Vector3d pixelFar(mWidth * sample.x(), mHeight * (1 - sample.y()), mFar);
        Eigen::Vector3d worldNear = unproject(pixelNear);
        Eigen::Vector3d worldFar  = unproject(pixelFar);
        Eigen::Vector3d dir       = worldFar - worldNear;
        Ray3d ray;
        ray.origin    = worldNear;
        ray.direction = dir.normalized();
        ray.tMax      = dir.norm();
        return { ray, Spectrum::Ones() };
    }

    std::pair<RayDifferential3d, Spectrum> ProjectiveCamera::sampleRayDifferential(const Eigen::Vector2d& sample) const
    {
        Eigen::Vector3d pixelNear(mWidth * sample.x(), mHeight * (1 - sample.y()), mNear);
        Eigen::Vector3d pixelFar(mWidth * sample.x(), mHeight * (1 - sample.y()), mFar);
        Eigen::Vector3d worldNear = unproject(pixelNear);
        Eigen::Vector3d worldFar  = unproject(pixelFar);
        Eigen::Vector3d dir       = worldFar - worldNear;

        RayDifferential3d ray(worldNear, dir.normalized());
        ray.tMax             = dir.norm();
        ray.hasDifferentials = true;

        Eigen::Vector3d worldNear_x = unproject(pixelNear + Eigen::Vector3d(1, 0, 0));
        Eigen::Vector3d worldFar_x  = unproject(pixelFar + Eigen::Vector3d(1, 0, 0));
        Eigen::Vector3d dir_x       = worldFar - worldNear;
        ray.origin_x                = worldNear_x;
        ray.direction_x             = (worldFar_x - worldNear_x).normalized();

        Eigen::Vector3d worldNear_y = unproject(pixelNear + Eigen::Vector3d(0, 1, 0));
        Eigen::Vector3d worldFar_y  = unproject(pixelFar + Eigen::Vector3d(0, 1, 0));
        Eigen::Vector3d dir_y       = worldFar - worldNear;
        ray.origin_y                = worldNear_y;
        ray.direction_y             = (worldFar_y - worldNear_y).normalized();

        return { ray, Spectrum::Ones() };
    }

    std::pair<DirectionSample, Spectrum> ProjectiveCamera::sampleDirection(const Interaction& it, const Eigen::Vector2d& sample) const
    {
        // this function will become relevant for bidirectional methods.
        throw std::logic_error("Not implemented yet.");
        return { DirectionSample(), Spectrum::Zero() };
    }
}
