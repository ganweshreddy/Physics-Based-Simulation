#pragma once

#include "endpoint.hpp"
#include "ray_differential.hpp"

#include <Eigen/Eigen>

namespace vislab
{
    /**
     * @brief Class that implements a basic camera.
     */
    class Camera : public Endpoint
    {
        VISLAB_INTERFACE(Camera, Endpoint)

    public:
        /**
         * @brief Constructor.
         */
        explicit Camera();

        /**
         * @brief Gets the view matrix.
         * @return View matrix.
         */
        const Eigen::Matrix4d& getView() const;

        /**
         * @brief Gets the position of the camera.
         * @return Camera position.
         */
        const Eigen::Vector3d& getPosition() const;

        /**
         * @brief Gets the location that the camera looks at.
         * @return Look at location.
         */
        const Eigen::Vector3d& getLookAt() const;

        /**
         * @brief Gets the up vector of the camera.
         * @return Up vector.
         */
        const Eigen::Vector3d& getUp() const;

        /**
         * @brief Sets the position of the camera.
         * @param Position Camera position.
         */
        void setPosition(const Eigen::Vector3d& Position);

        /**
         * @brief Sets the location that the camera looks at.
         * @param LookAt Look at location.
         */
        void setLookAt(const Eigen::Vector3d& LookAt);

        /**
         * @brief Sets the up vector of the camera.
         * @param Up Up vector.
         */
        void setUp(const Eigen::Vector3d& Up);

        /**
         * @brief Transforms all corners of a world space bounding box to view space and returns the view space bounding box of the transformed corners.
         * @param box World space bounding box to transform.
         * @return View space bounding box.
         */
        Eigen::AlignedBox3d getViewspaceAABB(const Eigen::AlignedBox3d& box) const;

        /**
         * @brief Samples a ray differential.
         * @param sample Sample in [0,1]^2.
         * @return Ray differential and importance.
         */
        virtual std::pair<RayDifferential3d, Spectrum> sampleRayDifferential(const Eigen::Vector2d& sample) const = 0;

    protected:
        /**
         * @brief Recomputes the view matrix.
         */
        void updateViewMatrix();

        /**
         * @brief View matrix.
         */
        Eigen::Matrix4d mView;

        /**
         * @brief Position of the camera.
         */
        Eigen::Vector3d mPosition;

        /**
         * @brief Location that the camera looks at.
         */
        Eigen::Vector3d mLookAt;

        /**
         * @brief Up vector of the camera.
         */
        Eigen::Vector3d mUp;
    };
}
