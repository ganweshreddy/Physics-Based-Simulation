#include <vislab/core/itransfer_function.hpp>

namespace vislab
{
    ITransferFunction::ITransferFunction(const double& minValue, const double& maxValue)
        : minValue(minValue)
        , maxValue(maxValue)
    {
    }
}
