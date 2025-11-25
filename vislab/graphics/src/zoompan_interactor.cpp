#include <vislab/graphics/zoompan_interactor.hpp>

#include <vislab/graphics/key_state.hpp>
#include <vislab/graphics/mouse_state.hpp>
#include <vislab/graphics/projective_camera.hpp>

namespace vislab
{
    ZoomPanInteractor::ZoomPanInteractor()
        : IInteractor()
        , mLastX(0)
        , mLastY(0)
        , mMovementSpeed(1)
    {
    }

    void ZoomPanInteractor::onMouseEvent(const MouseState& mouseState)
    {
        auto projectiveCamera = std::dynamic_pointer_cast<ProjectiveCamera>(mCamera);
        if (projectiveCamera == nullptr)
            return;

        int lastX = mLastX;
        int lastY = mLastY;
        if (mouseState.scrollDelta != 0)
        {
            double y = mouseState.scrollDelta;
            projectiveCamera->setWidth(projectiveCamera->getWidth() - y);
        }
        if (mouseState.middleIsDown && !mouseState.middleDown)
        {
            double scalefactor = projectiveCamera->getWidth();
            double mx          = (mouseState.x - lastX) * mMovementSpeed * scalefactor;
            double my          = (mouseState.y - lastY) * mMovementSpeed * scalefactor;
            if (mx != 0 || my != 0)
            {
                Eigen::Vector3d lookAt = projectiveCamera->getLookAt();
                Eigen::Vector3d eye    = projectiveCamera->getPosition();
                Eigen::Vector3d up     = projectiveCamera->getUp();
                Eigen::Vector3d dir    = lookAt - eye;
                Eigen::Vector3d right  = dir.cross(up).normalized();
                Eigen::Vector3d inup   = right.cross(dir).normalized();
                lookAt += right * mx + inup * my;
                eye += right * mx + inup * my;
                projectiveCamera->setLookAt(lookAt);
                projectiveCamera->setPosition(eye);
            }
        }
        mLastX = mouseState.x;
        mLastY = mouseState.y;
        mouseEvent.notify(this, &mouseState);
    }

    void ZoomPanInteractor::onKeyEvent(const KeyState& keyState)
    {
        keyEvent.notify(this, &keyState);
    }

    const double& ZoomPanInteractor::getMovementSpeed() const { return mMovementSpeed; }
    void ZoomPanInteractor::setMovementSpeed(const double& step) { mMovementSpeed = step; }
}
