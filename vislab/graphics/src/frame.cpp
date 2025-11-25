#include <vislab/graphics/frame.hpp>

namespace vislab
{
    Frame::Frame()
        : n(0, 0, 0)
        , s(0, 0, 0)
        , t(0, 0, 0)
    {
    }

    Frame::Frame(const Eigen::Vector3d& v)
        : n(v)
    {
        std::tie(s, t) = coordinate_system(v);
    }

    Frame::Frame(const Eigen::Vector3d& s, const Eigen::Vector3d& t, const Eigen::Vector3d& n)
        : s(s)
        , t(t)
        , n(n)
    {
    }

    Eigen::Vector3d Frame::toLocal(const Eigen::Vector3d& v) const
    {
        return Eigen::Vector3d(v.dot(s), v.dot(t), v.dot(n));
    }

    Eigen::Vector3d Frame::toWorld(const Eigen::Vector3d& v) const
    {
        return s * v.x() + t * v.y() + n * v.z();
    }

    double Frame::cos_theta(const Eigen::Vector3d& v) { return v.z(); }

    double Frame::cos_theta_2(const Eigen::Vector3d& v) { return (v.z() * v.z()); }

    double Frame::sin_theta(const Eigen::Vector3d& v) { return std::sqrt(std::max(0., sin_theta_2(v))); }

    double Frame::sin_theta_2(const Eigen::Vector3d& v) { return v.x() * v.x() + v.y() * v.y(); }

    double Frame::tan_theta(const Eigen::Vector3d& v)
    {
        double temp = 1. - v.z() * v.z();
        return std::sqrt(std::max(0., temp)) / v.z();
    }

    double Frame::tan_theta_2(const Eigen::Vector3d& v)
    {
        double temp = 1. - v.z() * v.z();
        return std::max(temp, 0.) / (v.z() * v.z());
    }

    double Frame::sin_phi(const Eigen::Vector3d& v)
    {
        double sin_theta_2   = Frame::sin_theta_2(v),
               inv_sin_theta = 1.0 / std::sqrt(Frame::sin_theta_2(v));
        return std::abs(sin_theta_2) <= (4. * std::numeric_limits<double>::epsilon()) ? 0 : std::min(std::max(-.1, v.y() * inv_sin_theta), 1.);
    }

    double Frame::cos_phi(const Eigen::Vector3d& v)
    {
        double sin_theta_2   = Frame::sin_theta_2(v),
               inv_sin_theta = 1.0 / std::sqrt(Frame::sin_theta_2(v));
        return std::abs(sin_theta_2) <= (4. * std::numeric_limits<double>::epsilon()) ? 0 : std::min(std::max(-.1, v.x() * inv_sin_theta), 1.);
    }

    std::pair<double, double> Frame::sincos_phi(const Eigen::Vector3d& v)
    {
        double sin_theta_2     = Frame::sin_theta_2(v),
               inv_sin_theta   = 1.0 / std::sqrt(Frame::sin_theta_2(v));
        Eigen::Vector2d result = Eigen::Vector2d(v.x(), v.y()) * inv_sin_theta;
        result                 = std::abs(sin_theta_2) <= (4. * std::numeric_limits<double>::epsilon()) ? Eigen::Vector2d(1, 0) : Eigen::Vector2d(std::min(std::max(-.1, result.x()), 1.), std::min(std::max(-.1, result.y()), 1.));
        return { result.y(), result.x() };
    }

    double Frame::sin_phi_2(const Eigen::Vector3d& v)
    {
        double sin_theta_2 = Frame::sin_theta_2(v);
        return std::abs(sin_theta_2) <= (4. * std::numeric_limits<double>::epsilon()) ? 0 : std::min(std::max(-.1, v.y() * v.y() / sin_theta_2), 1.);
    }

    double Frame::cos_phi_2(const Eigen::Vector3d& v)
    {
        double sin_theta_2 = Frame::sin_theta_2(v);
        return std::abs(sin_theta_2) <= (4. * std::numeric_limits<double>::epsilon()) ? 0 : std::min(std::max(-.1, v.x() * v.x() / sin_theta_2), 1.);
    }

    std::pair<double, double> Frame::sincos_phi_2(const Eigen::Vector3d& v)
    {
        double sin_theta_2     = Frame::sin_theta_2(v),
               inv_sin_theta_2 = 1.0 / (sin_theta_2);
        Eigen::Vector2d result = Eigen::Vector2d(v.x() * v.x(), v.y() * v.y()) * inv_sin_theta_2;
        result                 = std::abs(sin_theta_2) <= (4. * std::numeric_limits<double>::epsilon()) ? Eigen::Vector2d(1, 0) : Eigen::Vector2d(std::min(std::max(-.1, result.x()), 1.), std::min(std::max(-.1, result.y()), 1.));
        return { result.y(), result.x() };
    }

    bool Frame::operator==(const Frame& frame) const
    {
        return (frame.s == s) && (frame.t == t) && (frame.n == n);
    }

    bool Frame::operator!=(const Frame& frame) const
    {
        return (frame.s != s) || (frame.t != t) || (frame.n != n);
    }

    std::pair<Eigen::Vector3d, Eigen::Vector3d> Frame::coordinate_system(const Eigen::Vector3d& n)
    {
        double sign = n.z() < 0 ? -1 : 1,
               a    = -1.0 / (sign + n.z()),
               b    = n.x() * n.y() * a;
        return {
            Eigen::Vector3d(((n.x() * n.x()) * a * (n.z() < 0 ? -1 : 1)) + 1.f,
                            b * (n.z() < 0 ? -1 : 1),
                            n.x() * -(n.z() < 0 ? -1 : 1)),
            Eigen::Vector3d(b, sign + (n.y() * n.y()) * a, -n.y())
        };
    }
}
