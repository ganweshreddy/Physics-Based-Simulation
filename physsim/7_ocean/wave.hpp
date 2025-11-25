#include "dispersion.hpp"

#include <Eigen/Eigen>

namespace physsim
{
    /**
     * @brief Model of a single wave.
     */
    struct Wave
    {
        /**
         * @brief Wavelength lambda [m]
         */
        double waveLength;

        /**
         * @brief Amplitude A [m]
         */
        double amplitude;

        /**
         * @brief direction expressed by an angle [rad]
         */
        double angle;

        /**
         * @brief Steepness Q (Gerstner wave parameter) in [0,1].
         */
        double steepness;

        /**
         * @brief Phase phi
         */
        double phase;

        /**
         * @brief Dispersion model.
         */
        Dispersion dispersion;

        /**
         * @brief Wave number k. [rad/m]
         * @return k.
         */
        double waveNumber() const;

        /**
         * @brief Normalized wave traveling direction.
         * @return Normalized direction.
         */
        Eigen::Vector2d direction() const;

        /**
         * @brief Wave traveling direction scaled by the wave number k.
         * @return Wave vector kk.
         */
        Eigen::Vector2d waveVector() const;

        /**
         * @brief Computes the angular frequency omega = k * c_p [rad/s]
         * @return Angular frequency omega.
         */
        double angularFrequency() const;

        /**
         * @brief Evaluates the offset for this wave from the given position.
         * @param position Location at which to compute the wave offset.
         * @param t Time in seconds at which the wave is evaluated.
         * @return Offset due to this wave.
         */
        Eigen::Vector3d offset(const Eigen::Vector2d& position, const double& t) const;
    };
}
