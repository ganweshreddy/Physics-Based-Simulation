#include <vislab/graphics/math.hpp>

namespace vislab
{
    std::tuple<bool, double, double> solveQuadratic(const double& a, const double& b, const double& c)
    {
        // Is this perhaps a linear equation?
        bool linear_case = a == 0;

        // If so, we require b != 0
        bool valid_linear = linear_case && b != 0;

        // Initialize solution with that of linear equation
        double x0, x1;
        x0 = x1 = -c / b;

        // Check if the quadratic equation is solvable
        double discrim       = b * b - 4 * a * c;
        bool valid_quadratic = !linear_case && (discrim >= 0);

        if (valid_quadratic)
        {
            double sqrt_discrim = std::sqrt(discrim);

            /* Numerically stable version of (-b (+/-) sqrt_discrim) / (2 * a)
             *
             * Based on the observation that one solution is always
             * accurate while the other is not. Finds the solution of
             * greater magnitude which does not suffer from loss of
             * precision and then uses the identity x1 * x2 = c / a
             */
            double temp = -0.5 * (b + std::copysign(sqrt_discrim, b));

            double x0p = temp / a,
                   x1p = c / temp;

            // Order the results so that x0 < x1
            double x0m = std::min(x0p, x1p),
                   x1m = std::max(x0p, x1p);

            x0 = linear_case ? x0 : x0m;
            x1 = linear_case ? x0 : x1m;
        }

        return std::make_tuple(valid_linear || valid_quadratic, x0, x1);
    }

    double unitAngleZ(const Eigen::Vector3d& v)
    {
        double temp = 2.f * std::asin(.5 * std::sqrt(v.x() * v.x() + v.y() * v.y() +
                                                     (v.z() - copysign(1., v.z())) * (v.z() - copysign(1., v.z()))));
        return v.z() >= 0 ? temp : (EIGEN_PI)-temp;
    }

    double degree2Radians(double d)
    {
        return d * (EIGEN_PI / 180.0);
    }

    double radians2Degree(double r)
    {
        return r * (180.0 / EIGEN_PI);
    }

    double safe_sqrt(double x) { return std::sqrt(std::max(0., x)); }

    double epsilon() { return std::numeric_limits<double>::epsilon() / 2; }

    double rayEpsilon() { return epsilon() * 1500; }

    double shadowEpsilon() { return rayEpsilon() * 10; }

    std::pair<Eigen::Vector3d, Eigen::Vector3d> coordinateSystem(const Eigen::Vector3d& n)
    {
        using Float = Eigen::Vector3d::Scalar;

        /* Based on "Building an Orthonormal Basis, Revisited" by
           Tom Duff, James Burgess, Per Christensen,
           Christophe Hery, Andrew Kensler, Max Liani,
           and Ryusuke Villemin (JCGT Vol 6, No 1, 2017) */

        Float sign = n.z() < 0 ? -1 : 1,
              a    = -1. / (sign + n.z()),
              b    = n.x() * n.y() * a;

        return {
            Eigen::Vector3d(((n.x() * n.x()) * a * (n.z() < 0 ? -1 : 1)) + 1.,
                            (b * (n.z() < 0 ? -1 : 1)),
                            (n.x() * -(n.z() < 0 ? -1 : 1))),
            Eigen::Vector3d(b, sign + n.y() * (n.y()) * a, -n.y())
        };
    }
}
