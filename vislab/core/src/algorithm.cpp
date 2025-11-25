#include <vislab/core/algorithm.hpp>

#include <stdexcept>

namespace vislab
{
    Algorithm::Algorithm()
    {
        active.setValue(true, false);
    }

    Algorithm::~Algorithm()
    {
    }

    UpdateInfo Algorithm::update(ProgressInfo& progress)
    {
        return internalUpdate(progress);
    }

    UpdateInfo Algorithm::update()
    {
        ProgressInfo progressInfo;
        return update(progressInfo);
    }
}
