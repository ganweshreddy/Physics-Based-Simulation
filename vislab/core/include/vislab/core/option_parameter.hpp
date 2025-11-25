#pragma once

#include "parameter.hpp"

namespace vislab
{
    /**
     * @brief Parameter that selects an option from a list of values (strings).
     */
    class OptionParameter : public Parameter<int32_t>
    {
        VISLAB_OBJECT(OptionParameter, Parameter<int32_t>)

    public:
        /**
         * @brief Constructor.
         */
        OptionParameter();

        /**
         * @brief Gets the list of labels.
         * @return Vector with labels for the different options.
         */
        [[nodiscard]] const std::vector<std::string>& getLabels() const;

        /**
         * @brief Sets the list of labels.
         * @param newLabels New list of labels to set.
         * @param notifyOnChange Flag that determines whether an event is raised to all listeners.
         */
        void setLabels(const std::vector<std::string>& newLabels, bool notifyOnChange = true);

        /**
         * @brief Gets the currently selected label.
         * @return The currently selected label.
         */
        std::string getCurrentLabel() const;

        /**
         * @brief Tests if another parameter stores the exact same values and whether they have the same visibilty.
         * @param other Other parameter to compare with.
         * @return True, if the values and visibility are the same.
         */
        bool operator==(const OptionParameter& other) const;

        /**
         * @brief Event to be raised when the labels change.
         */
        StringArrayEvent onLabelsChange;

    private:
        /**
         * @brief Vector containing the labels for the different options.
         */
        std::vector<std::string> mLabels;
    };
}
