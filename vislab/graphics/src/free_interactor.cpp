#include <vislab/graphics/camera.hpp>
#include <vislab/graphics/free_interactor.hpp>
#include <vislab/graphics/key_state.hpp>
#include <vislab/graphics/mouse_state.hpp>

namespace vislab
{
    FreeInteractor::FreeInteractor()
        : IInteractor()
        , mLastX(0)
        , mLastY(0)
        , mMovementSpeed(0.1f)
        , mRotationSpeed(0.001f)
    {
    }

    void FreeInteractor::onMouseEvent(const MouseState& mouseState)
    {
        if (mCamera == nullptr)
            return;
        int lastX = mLastX;
        int lastY = mLastY;
        if (mouseState.rightIsDown)
        {
            double rx = (mouseState.x - lastX) * mRotationSpeed;
            double ry = -(mouseState.y - lastY) * mRotationSpeed;
            if (rx != 0 || ry != 0)
            {
                Eigen::Vector3d lookAt = mCamera->getLookAt();
                Eigen::Vector3d eye    = mCamera->getPosition();
                Eigen::Vector3d up     = mCamera->getUp();

                // rotate left - right
                Eigen::Vector3d dir = lookAt - eye;
                dir.normalize();
                Eigen::Quaterniond rot = Eigen::Quaterniond(Eigen::AngleAxisd(rx, up));
                dir                    = rot * dir;

                // rotate up - down
                Eigen::Vector3d lastDirection = dir;
                Eigen::Vector3d right         = dir.cross(up).normalized();
                rot                           = Eigen::Quaterniond(Eigen::AngleAxisd(ry, right));
                dir                           = rot * dir;

                // prevent flipping around the up-vector
                if (std::abs(Eigen::Vector3d::dot(dir, up)) > 0.99)
                    dir = lastDirection;

                lookAt = eye + dir;
                mCamera->setLookAt(lookAt);
                mCamera->setPosition(eye);
            }
        }
        mLastX = mouseState.x;
        mLastY = mouseState.y;
        mouseEvent.notify(this, &mouseState);
    }

    void FreeInteractor::onKeyEvent(const KeyState& keyState)
    {
        if (mCamera == nullptr)
            return;
        double x = ((keyState.isDown_A ? 1 : 0) - (keyState.isDown_D ? 1 : 0)) * mMovementSpeed;
        double y = ((keyState.isDown_W ? 1 : 0) - (keyState.isDown_S ? 1 : 0)) * mMovementSpeed;
        if (x != 0 || y != 0)
        {
            Eigen::Vector3d lookAt = mCamera->getLookAt();
            Eigen::Vector3d eye    = mCamera->getPosition();
            Eigen::Vector3d up     = mCamera->getUp();
            Eigen::Vector3d dir    = lookAt - eye;
            if (y != 0)
            {
                eye    = eye + dir * y;
                lookAt = lookAt + dir * y;
            }
            if (x != 0)
            {
                Eigen::Vector3d right = dir.cross(up).normalized();
                eye                   = eye + right * x;
                lookAt                = lookAt + right * x;
            }
            mCamera->setLookAt(lookAt);
            mCamera->setPosition(eye);
        }
        keyEvent.notify(this, &keyState);
    }

    const double& FreeInteractor::getMovementSpeed() const { return mMovementSpeed; }
    void FreeInteractor::setMovementSpeed(const double& step) { mMovementSpeed = step; }
    const double& FreeInteractor::getRotationSpeed() const { return mRotationSpeed; }
    void FreeInteractor::setRotationSpeed(const double& step) { mRotationSpeed = step; }
}
