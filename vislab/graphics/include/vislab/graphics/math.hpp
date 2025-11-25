#pragma once

#include <Eigen/Eigen>

namespace vislab
{
    /**
     * @brief Solve a quadratic equation of the form a* x ^ 2 + b* x + c = 0. his function is adapted from https://github.com/mitsuba-renderer/mitsuba2
     * @param a Coefficient of the quadratic monomial.
     * @param b Coefficient of the linear monomial.
     * @param c Coefficient of the constant part.
     * @return Tuple containing true if there are real roots, as well as the roots themselves sorted in ascending order.
     */
    std::tuple<bool, double, double> solveQuadratic(const double& a, const double& b, const double& c);

    /**
     * @brief Numerically well-behaved routine for computing the angle between the unit direction vector 'v' and the z-axis. This should be used wherever one is tempted to compute std::acos(v.z()). By Don Hatch at http://www.plunk.org/~hatch/rightway.php
     * @param v Unit vector.
     * @return Enclosed angle.
     */
    double unitAngleZ(const Eigen::Vector3d& v);

    /**
     * @brief Converts an angle from degree to radians.
     * @param d Angle in degree.
     * @return Angle in radians.
     */
    double degree2Radians(double d);

    /**
     * @brief Converts an angle from radiants to degree.
     * @param r Angle in radians.
     * @return Angle in degree.
     */
    double radians2Degree(double r);

    /**
     * @brief Square root of number. Clamped to zero if number is negative.
     * @param x Number to take square root of.
     * @return Square root.
     */
    double safe_sqrt(double x);

    /**
     * @brief Numerical epsilon.
     * @return Epsilon for numbers close to zero.
     */
    double epsilon();

    /**
     * @brief Global ray offset epsilon.
     * @return Ray offset for avoidance of self intersections.
     */
    double rayEpsilon();

    /**
     * @brief Global shadow ray offset epsilon.
     * @return Ray offset for avoidance of self intersections.
     */
    double shadowEpsilon();

    /**
     * @brief Complete the set {a} to an orthonormal basis {a, b, c}
     * @param n Normal to constract basis around.
     * @return The other two basis vectors.
     */
    std::pair<Eigen::Vector3d, Eigen::Vector3d> coordinateSystem(const Eigen::Vector3d& n);
}
