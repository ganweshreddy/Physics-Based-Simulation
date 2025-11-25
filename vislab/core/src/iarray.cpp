#include <vislab/core/iarray.hpp>

namespace vislab
{
    IArray::IArray()
        : mName("")
    {
    }

    const std::string& IArray::getName() const
    {
        return mName;
    }

    void IArray::setName(const std::string& name)
    {
        mName = name;
    }

    Eigen::VectorXd IArray::getValueDouble(const Eigen::Index& index) const
    {
        return getValueDoubleImpl(index);
    }

    void IArray::setValueDouble(const Eigen::Index& index, const Eigen::VectorXd& value)
    {
        setValueDoubleImpl(index, value);
    }

    Eigen::VectorXl IArray::getValueInt64(const Eigen::Index& index) const
    {
        return getValueInt64Impl(index);
    }

    void IArray::setValueInt64(const Eigen::Index& index, const Eigen::VectorXl& value)
    {
        setValueInt64Impl(index, value);
    }
}
