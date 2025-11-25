#include "wave.hpp"

namespace physsim
{
    double Wave::waveNumber() const { return 2 * (float)EIGEN_PI / waveLength; }

    Eigen::Vector2d Wave::direction() const { return Eigen::Vector2d(std::cos(angle), std::sin(angle)); }

    Eigen::Vector2d Wave::waveVector() const { return direction() * waveNumber(); }

    double Wave::angularFrequency() const
    {
        double k = waveNumber();
        return k * dispersion.phaseSpeed(k);
    }

    Eigen::Vector3d Wave::offset(const Eigen::Vector2d& position, const double& t) const
    {
        // TODO: compute the offset for the displaced wave from the given position.
        return Eigen::Vector3d(0, 0, 0);
    }
}
