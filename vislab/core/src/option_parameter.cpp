#include <vislab/core/option_parameter.hpp>

namespace vislab
{
    OptionParameter::OptionParameter()
        : Parameter<int32_t>()
    {
    }

    const std::vector<std::string>& OptionParameter::getLabels() const
    {
        return mLabels;
    }

    void OptionParameter::setLabels(const std::vector<std::string>& newLabels, bool notifyOnChange)
    {
        mLabels = newLabels;
        if (notifyOnChange)
            onLabelsChange.notify(this, &mLabels);
    }

    std::string OptionParameter::getCurrentLabel() const
    {
        return mLabels.at(mValue);
    }

    bool OptionParameter::operator==(const OptionParameter& other) const
    {
        return Parameter<int32_t>::operator==(other) && mLabels == other.mLabels;
    }
}
