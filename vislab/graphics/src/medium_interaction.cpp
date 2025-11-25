#include <vislab/graphics/medium_interaction.hpp>

namespace vislab
{
    Eigen::Vector3d MediumInteraction::toWorld(const Eigen::Vector3d& v) const
    {
        return sh_frame.toWorld(v);
    }

    Eigen::Vector3d MediumInteraction::toLocal(const Eigen::Vector3d& v) const
    {
        return sh_frame.toLocal(v);
    }
}
