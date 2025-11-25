#pragma once

#include "camera.hpp"

namespace vislab
{
    /**
     * @brief Class that implements a projective camera.
     */
    class ProjectiveCamera : public Camera
    {
        VISLAB_INTERFACE(ProjectiveCamera, Camera)

    public:
        /**
         * @brief Constructor.
         */
        explicit ProjectiveCamera();

        /**
         * @brief Gets the projection matrix.
         * @return Projection matrix.
         */
        const Eigen::Matrix4d& getProj() const;

        /**
         * @brief Gets the near plane distance of the viewing frustum.
         * @return Near plane distance.
         */
        const double& getNear() const;

        /**
         * @brief Gets the far plane distance of the viewing frustum.
         * @return Far plane distance.
         */
        const double& getFar() const;

        /**
         * @brief Gets the width of the viewport.
         * @return Width of the viewport.
         */
        const double& getWidth() const;

        /**
         * @brief Gets the height of the viewport.
         * @return Height of the viewport.
         */
        const double& getHeight() const;

        /**
         * @brief Gets the aspect ratio of the projection.
         * @return Aspect ratio of the viewport.
         */
        double getAspectRatio() const;

        /**
         * @brief Sets the near plane distance of the viewing frustum.
         * @param  Near plane distance.
         */
        void setNear(const double& near);

        /**
         * @brief Sets the far plane distance of the viewing frustum.
         * @param  Far plane distance.
         */
        void setFar(const double& far);

        /**
         * @brief Sets the width of the viewport.
         * @param width Width of the viewport.
         */
        void setWidth(const double& width);

        /**
         * @brief Sets the height of the viewport.
         * @param height Height of the viewport.
         */
        void setHeight(const double& height);

        /**
         * @brief Projects a 3D world location into screen space coordinates.
         * @param coord World coordinate to project into screen space.
         * @return Screen space coordinate.
         */
        Eigen::Vector2d project(const Eigen::Vector3d& coord) const;

        /**
         * @brief Projects a screen space coordinate (with depth) back into world space.
         * @param coord Screen space coordinate to unproject into world space.
         * @return World space coordinate.
         */
        Eigen::Vector3d unproject(const Eigen::Vector3d& coord) const;

        /**
         * @brief Generates a ray in world space.
         * @param sample Sample location in [0,1]^2.
         * @param sample2 Another sample location in [0,1]^2.
         * @return Ray that was sampled and spectrum.
         */
        std::pair<Ray3d, Spectrum> sampleRay(const Eigen::Vector2d& sample, const Eigen::Vector2d& sample2) const override;

        /**
         * @brief Samples a ray differential.
         * @param sample Sample in [0,1]^2.
         * @return Ray differential and importance.
         */
        std::pair<RayDifferential3d, Spectrum> sampleRayDifferential(const Eigen::Vector2d& sample) const override;

        /**
         * @brief Given a reference point in the scene, sample a direction from the reference point towards the endpoint (ideally proportional to the emission/sensitivity profile).
         * @param it Interaction to sample the direction from.
         * @param sample Sample location in [0,1]^2.
         * @return Information about the sampled direction.
         */
        std::pair<DirectionSample, Spectrum> sampleDirection(const Interaction& it, const Eigen::Vector2d& sample) const override;

    protected:
        /**
         * @brief Recomputes the projection matrix.
         */
        virtual void updateProjMatrix() = 0;

        /**
         * @brief Near plane distance of the viewing frustum.
         */
        double mNear;

        /**
         * @brief Far plane distance of the viewing frustum.
         */
        double mFar;

        /**
         * @brief Width of the viewport.
         */
        double mWidth;

        /**
         * @brief Height of the viewport.
         */
        double mHeight;

        /**
         * @brief Projection matrix.
         */
        Eigen::Matrix4d mProj;
    };
}
