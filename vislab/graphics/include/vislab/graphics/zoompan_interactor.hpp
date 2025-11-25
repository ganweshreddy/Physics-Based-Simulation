#pragma once

#include "iinteractor.hpp"

namespace vislab
{
    /**
     * @brief Interactor that implements a zooming and panning camera.
     */
    class ZoomPanInteractor : public IInteractor
    {
        VISLAB_OBJECT(ZoomPanInteractor, IInteractor)

    public:
        /**
         * @brief Constructor.
         */
        explicit ZoomPanInteractor();

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

    private:
        /**
         * @brief Scaling factor for the movement speed.
         */
        double mMovementSpeed;

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
