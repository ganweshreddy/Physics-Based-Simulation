#pragma once

#include <Eigen/Eigen>

namespace vislab
{
    /**
     * @brief Class that represents a local coordinate system.
     */
    class Frame
    {
    public:
        /**
         * @brief Tangent vector for s-direction.
         */
        Eigen::Vector3d s;

        /**
         * @brief Bi-Tangent vector for t-direction.
         */
        Eigen::Vector3d t;

        /**
         * @brief Normal of the frame.
         */
        Eigen::Vector3d n;

        /**
         * @brief Empty constructor for zero frame.
         */
        Frame();

        /**
         * @brief Construct a new coordinate frame from a given normal.
         * @param v Normal to construct shading frame around.
         */
        Frame(const Eigen::Vector3d& v);

        /**
         * @brief Constructs a frame from three given basis vectors.
         * @param s s-tangent.
         * @param t t-tangent.
         * @param n Normal.
         */
        Frame(const Eigen::Vector3d& s, const Eigen::Vector3d& t, const Eigen::Vector3d& n);

        /**
         * @brief Convert from world coordinates to local coordinates.
         * @param v Vector to transform.
         * @return Transformed vector.
         */
        Eigen::Vector3d toLocal(const Eigen::Vector3d& v) const;

        /**
         * @brief Convert from local coordinates to world coordinates.
         * @param v Vector to transform.
         * @return Transformed vector.
         */
        Eigen::Vector3d toWorld(const Eigen::Vector3d& v) const;

        /**
         * @brief Given a unit direction, return the cosine of the elevation angle in a reference spherical coordinate system.
         * @param v Given unit vector.
         * @return Cosine of the elevation angle.
         */
        static double cos_theta(const Eigen::Vector3d& v);

        /**
         * @brief Given a unit direction, return the square cosine of the elevation angle in a reference spherical coordinate system.
         * @param v Given unit vector.
         * @return Square cosine of the elevation angle.
         */
        static double cos_theta_2(const Eigen::Vector3d& v);

        /**
         * @brief Given a unit direction, return the sine of the elevation angle in a reference spherical coordinate system.
         * @param v Given unit vector.
         * @return Sine of the elevation angle.
         */
        static double sin_theta(const Eigen::Vector3d& v);

        /**
         * @brief Given a unit direction, return the square sine of the elevation angle in a reference spherical coordinate system.
         * @param v Given unit vector.
         * @return Square sine of the elevation angle.
         */
        static double sin_theta_2(const Eigen::Vector3d& v);

        /**
         * @brief Given a unit direction, return the tangent of the elevation angle in a reference spherical coordinate system.
         * @param v Given unit vector.
         * @return Tangent of the elevation angle.
         */
        static double tan_theta(const Eigen::Vector3d& v);

        /**
         * @brief Given a unit direction, return the square tangent of the elevation angle in a reference spherical coordinate system.
         * @param v Given unit vector.
         * @return Square tangent of the elevation angle.
         */
        static double tan_theta_2(const Eigen::Vector3d& v);

        /**
         * @brief Give a unit direction, return the sine of the azimuth in a reference spherical coordinate system.
         * @param v Given unit vector.
         * @return Sine of the azimuth.
         */
        static double sin_phi(const Eigen::Vector3d& v);

        /**
         * @brief Given a unit direction, return the cosine of the azimuth in a reference spherical coordinate system.
         * @param v Given unit vector.
         * @return Cosine of the azimuth.
         */
        static double cos_phi(const Eigen::Vector3d& v);

        /**
         * @brief Given a unit direction, return the sine and cosine of the azimuth in a reference spherical coordinate system.
         * @param v Given unit vector.
         * @return Sine and cosine of the azimuth.
         */
        static std::pair<double, double> sincos_phi(const Eigen::Vector3d& v);

        /**
         * @brief Given a unit direction, return the squared sine of the azimuth in a reference spherical coordinate system.
         * @param v Given unit vector.
         * @return Squared sine of the azimuth.
         */
        static double sin_phi_2(const Eigen::Vector3d& v);

        /**
         * @brief Given a unit direction, return the squared cosine of the azimuth in a reference spherical coordinate system.
         * @param v Given unit vector.
         * @return Squared cosine of the azimuth.
         */
        static double cos_phi_2(const Eigen::Vector3d& v);

        /**
         * @brief Given a unit direction, return the squared sine and cosine of the azimuth in a reference spherical coordinate system.
         * @param v Given unit vector.
         * @return Squared sine and cosine of the azimuth.
         */
        static std::pair<double, double> sincos_phi_2(const Eigen::Vector3d& v);

        /**
         * @brief Equality test.
         * @param frame Frame to compare to.
         * @return True if equal.
         */
        bool operator==(const Frame& frame) const;

        /**
         * @brief Inequality test.
         * @param frame Frame to compare to.
         * @return True if unequal.
         */
        bool operator!=(const Frame& frame) const;

    private:
        /**
         * @brief Based on "Building an Orthonormal Basis, Revisited" by Tom Duff, James Burgess, Per Christensen, Christophe Hery, Andrew Kensler, Max Liani, and Ryusuke Villemin (JCGT Vol 6, No 1, 2017)
         * @param n Given normal vector.
         * @return Tangent and bi-tangent.
         */
        std::pair<Eigen::Vector3d, Eigen::Vector3d> coordinate_system(const Eigen::Vector3d& n);
    };
}
