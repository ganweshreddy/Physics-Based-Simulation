#pragma once

#include "data.hpp"
#include "event.hpp"
#include "iparameter.hpp"
#include "parameter_fwd.hpp"
#include "traits.hpp"

namespace vislab
{
    /**
     * @brief Base class for a parameter.
     * @tparam TValueType Type of value that is stored in the parameter.
     */
    template <typename TValueType>
    class Parameter : public IParameter
    {
        VISLAB_OBJECT(Parameter<TValueType>, IParameter)

    public:
        /**
         * @brief Type of the value stored in this parameter.
         */
        using ValueType = TValueType;

        /**
         * @brief Constructor.
         */
        Parameter()
            : mHidden(false)
        {
            if constexpr (std::is_scalar_v<TValueType>)
                mValue = TValueType(0);
            else if constexpr (detail::is_eigen_v<TValueType>)
                mValue = TValueType::Zero();
            else
                mValue = TValueType();
        }

        /**
         * @brief Constructor.
         * @param value Initial value.
         */
        Parameter(const TValueType& value)
            : mHidden(false)
            , mValue(value)
        {
        }

        /**
         * @brief Copy-assignment operator. Event listeners are not copied.
         * @param other Other parameter to copy from.
         * @return Reference to self.
         */
        Parameter& operator=(const Parameter& other)
        {
            mValue  = other.mValue;
            mHidden = other.mHidden;
            return *this;
        }

        /**
         * @brief Gets the value stored in this parameter.
         * @return Values stored in this parameter.
         */
        const TValueType& getValue() const { return mValue; }

        /**
         * @brief Sets the value stored in this parameter.
         * @param value New value to be set.
         * @param notifyOnChange Flag that determines whether an event is raised to all listeners.
         */
        virtual void setValue(const TValueType& value, bool notifyOnChange = true)
        {
            if (mValue != value)
            {
                mValue = value;
                if (notifyOnChange)
                    onChange.notify(this, &mValue);
            }
        }

        /**
         * @brief Event that is raised when the value changes.
         */
        TEvent<Parameter, TValueType> onChange;

        /**
         * @brief Is the parameter hidden from the UI?
         * @return True if the parameter is hidden from the UI.
         */
        bool isHidden() const { return mHidden; }

        /**
         * @brief Tests if another parameter stores the exact same value and whether they have the same visibilty.
         * @param other Other parameter to compare with.
         * @return True, if the values and visibility are the same.
         */
        bool operator==(const Parameter& other) const
        {
            return mValue == other.mValue && mHidden == other.mHidden;
        }

        /**
         * @brief Hides this parameter.
         */
        void hide()
        {
            if (!mHidden)
            {
                mHidden = true;
                onHiddenChange.notify(this, &mHidden);
            }
        }

        /**
         * @brief Shows this parameter.
         */
        void show()
        {
            if (mHidden)
            {
                mHidden = false;
                onHiddenChange.notify(this, &mHidden);
            }
        }

        /**
         * @brief Event that is raised when the hidden flag changes.
         */
        BoolEvent onHiddenChange;

    protected:
        /**
         * @brief Data value stored by this parameter.
         */
        TValueType mValue;

        /**
         * @brief Flag that stored whether this parameter is hidden from the UI.
         */
        bool mHidden;
    };
}
