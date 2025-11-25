#pragma once

#include "projective_camera.hpp"

namespace vislab
{
    /**
     * @brief Class that implements a projective camera.
     */
    class PerspectiveCamera : public ProjectiveCamera
    {
        VISLAB_OBJECT(PerspectiveCamera, ProjectiveCamera)

    public:
        /**
         * @brief Constructor.
         */
        explicit PerspectiveCamera();

        /**
         * @brief Gets the horizontal field of view (in radians) for perspective projection.
         * @return Horizontal field of view (in radians).
         */
        const double& getHorizontalFieldOfView() const;

        /**
         * @brief Sets the horizontal field of view (in radians) for perspective projection.
         * @param hfov Horizontal field of view (in radians).
         */
        void setHorizontalFieldOfView(const double& hfov);

    private:
        /**
         * @brief Recomputes the projection matrix.
         */
        void updateProjMatrix() override;

        /**
         * @brief Horizontal field of view (in radians) for perspective projection.
         */
        double mHorizontalFieldOfView;
    };
}
