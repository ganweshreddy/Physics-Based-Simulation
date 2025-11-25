#pragma once

#include "iinteractor.hpp"

namespace vislab
{
    /**
     * @brief Interactor that implements a trackball camera.
     */
    class TrackballInteractor : public IInteractor
    {
        VISLAB_OBJECT(TrackballInteractor, IInteractor)

    public:
        /**
         * @brief Constructor.
         */
        explicit TrackballInteractor();

        /**
         * @brief Processes mouse events.
         * @param mouseState Current mouse state.
         */
        void onMouseEvent(const MouseState& mouseState) override;

        /**
         * @brief Processes key events.
         * @param keyState Current key state.
         */
        void onKeyEvent(const KeyState& keyState) override;

        /**
         * @brief Gets the scaling factor for the movement speed.
         * @return Movement speed.
         */
        [[nodiscard]] const double& getMovementSpeed() const;

        /**
         * @brief Sets the scaling factor for the movement speed.
         * @param scale Movement speed.
         */
        void setMovementSpeed(const double& scale);

        /**
         * @brief Gets the scaling factor for the rotation speed.
         * @return Rotation speed.
         */
        [[nodiscard]] const double& getRotationSpeed() const;

        /**
         * @brief Sets the scaling factor for the rotation speed.
         * @param scale Rotation speed.
         */
        void setRotationSpeed(const double& scale);

    private:
        /**
         * @brief Scaling factor for the movement speed.
         */
        double mMovementSpeed;

        /**
         * @brief Scaling factor for the rotation speed.
         */
        double mRotationSpeed;

        /**
         * @brief Last x coordinate of the cursor.
         */
        int mLastX;

        /**
         * @brief Last y coordinate of the cursor.
         */
        int mLastY;
    };
}
