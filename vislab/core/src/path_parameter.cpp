#include <vislab/core/path_parameter.hpp>

namespace vislab
{
    PathParameter::PathParameter()
        : Parameter<std::string>()
        , mFilter("")
        , mFile(EFile::In)
    {
    }

    PathParameter::PathParameter(const EFile& fileDirection, const std::string& filter)
        : Parameter<std::string>()
        , mFilter(filter)
        , mFile(fileDirection)
    {
    }

    const PathParameter::EFile& PathParameter::getFileDirection() const
    {
        return mFile;
    }

    void PathParameter::setFileDirection(const PathParameter::EFile& file)
    {
        mFile = file;
    }

    const std::string& PathParameter::getFilter() const
    {
        return mFilter;
    }

    void PathParameter::setFilter(const std::string& filter)
    {
        mFilter = filter;
    }

    bool PathParameter::operator==(const PathParameter& other) const
    {
        return Parameter<std::string>::operator==(other) && mFilter == other.mFilter && mFile == other.mFile;
    }
}
