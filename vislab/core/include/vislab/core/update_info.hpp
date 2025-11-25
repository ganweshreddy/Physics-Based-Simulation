#pragma once

#include <string>

namespace vislab
{
    /**
     * @brief Return type for the Update function, which contains some general information about the success and a potential error message.
     */
    class UpdateInfo
    {
    public:
        /**
         * @brief State after the Update call.
         */
        enum class EValidationState
        {
            /**
             * @brief The update function terminated as expected.
             */
            Valid,

            /**
             * @brief A warning was reported but the update function was still able to terminate.
             */
            Warning,

            /**
             * @brief An error occured and the result of the update function is not valid.
             */
            Error
        };

        /**
         * @brief The Update call was successful.
         * @return Reports that a computation was valid.
         */
        static UpdateInfo reportValid();

        /**
         * @brief Report a warning message.
         * @param string Warning message.
         * @return Reports that a computation finished with a warning.
         */
        static UpdateInfo reportWarning(const std::string& string);

        /**
         * @brief Report an error message.
         * @param string Error message.
         * @return Reports that a computation failed due to an error.
         */
        static UpdateInfo reportError(const std::string& string);

        /**
         * @brief Gets the validation state.
         * @return Last validation state.
         */
        const EValidationState& getState() const;

        /**
         * @brief Gets the warning or error message.
         * @return Last error message.
         */
        const std::string& getError() const;

        /**
         * @brief Returns true if the computation was a success.
         * @return True if GetState() == UpdateInfo::EValidationState::Valid.
         */
        bool success() const;

    private:
        /**
         * @brief Constructor for an update info.
         * @param state Validation state.
         * @param string Optional error message.
         */
        UpdateInfo(const EValidationState& state, const std::string& string);

        /**
         * @brief State after the Update call.
         */
        EValidationState mValidationState;

        /**
         * @brief Optional warning or error message.
         */
        std::string mValidationError;
    };
}
