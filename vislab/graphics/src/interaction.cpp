#include <vislab/graphics/interaction.hpp>

namespace vislab
{
    // ======================================================================
    Interaction::Interaction()
        : t(0)
        , position(0, 0, 0)
    {
    }

    Interaction::Interaction(double tt, const Eigen::Vector3d& pos)
        : t(tt)
        , position(pos)
    {
    }

    bool Interaction::isValid() const
    {
        return t != std::numeric_limits<double>::infinity();
    }
}
