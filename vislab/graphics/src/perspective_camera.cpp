#include <vislab/graphics/perspective_camera.hpp>

namespace vislab
{
    PerspectiveCamera::PerspectiveCamera()
        : mHorizontalFieldOfView(EIGEN_PI / 4)
    {
        updateProjMatrix();
    }

    const double& PerspectiveCamera::getHorizontalFieldOfView() const { return mHorizontalFieldOfView; }

    void PerspectiveCamera::setHorizontalFieldOfView(const double& hfov)
    {
        mHorizontalFieldOfView = hfov;
        updateProjMatrix();
    }

    void PerspectiveCamera::updateProjMatrix()
    {
        mProj = Eigen::Matrix4d::perspectiveFovLH(mHorizontalFieldOfView, getAspectRatio(), mNear, mFar);
    }
}
