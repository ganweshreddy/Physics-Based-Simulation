#include "vislab/graphics/frame.hpp"

#include "Eigen/Eigen"
#include "gtest/gtest.h"

namespace vislab
{
    TEST(graphics, frame)
    {
        double thetas[2] = {
            30 * EIGEN_PI / 180,
            95 * EIGEN_PI / 180
        };

        for (int i = 0; i < 2; ++i)
        {
            double theta = thetas[i];
            double phi = 73 * EIGEN_PI / 180;
            double sin_theta = std::sin(theta);
            double cos_theta = std::cos(theta);
            double sin_phi = std::sin(phi);
            double cos_phi = std::cos(phi);
            Eigen::Vector3d v(
                cos_phi * sin_theta,
                sin_phi * sin_theta,
                cos_theta);
            Frame f(Eigen::Vector3d(1.0, 2.0, 3.0) / std::sqrt(14));
            Eigen::Vector3d v2 = f.toLocal(v);
            Eigen::Vector3d v3 = f.toWorld(v2);

            EXPECT_NEAR((v3 - v).squaredNorm(), 0, 1E-7);
            EXPECT_NEAR(Frame::cos_theta(v), cos_theta, 1E-7);
            EXPECT_NEAR(Frame::sin_theta(v), sin_theta, 1E-7);
            EXPECT_NEAR(Frame::cos_phi(v), cos_phi, 1E-7);
            EXPECT_NEAR(Frame::sin_phi(v), sin_phi, 1E-7);
            EXPECT_NEAR(Frame::cos_theta_2(v), cos_theta * cos_theta, 1E-7);
            EXPECT_NEAR(Frame::sin_theta_2(v), sin_theta * sin_theta, 1E-7);
            EXPECT_NEAR(Frame::cos_phi_2(v), cos_phi * cos_phi, 1E-7);
            EXPECT_NEAR(Frame::sin_phi_2(v), sin_phi * sin_phi, 1E-7);
            EXPECT_NEAR(Frame::sincos_phi(v).first, sin_phi, 1E-7);
            EXPECT_NEAR(Frame::sincos_phi(v).second, cos_phi, 1E-7);
            EXPECT_NEAR(Frame::sincos_phi_2(v).first, sin_phi * sin_phi, 1E-7);
            EXPECT_NEAR(Frame::sincos_phi_2(v).second, cos_phi * cos_phi, 1E-7);
        }
    }
}
