#include "vislab/graphics/fresnel.hpp"

#include "Eigen/Eigen"
#include "gtest/gtest.h"

namespace vislab
{
    TEST(graphics, fresnel)
    {
        using namespace vislab;

        double ct_crit = -sqrt(1 - 1. / 2.25);

        EXPECT_NEAR(std::get<0>(fresnel(1, 1.5)), 0.04, 1E-7);
        EXPECT_NEAR(std::get<1>(fresnel(1, 1.5)), -1, 1E-7);
        EXPECT_NEAR(std::get<2>(fresnel(1, 1.5)), 1.5, 1E-7);
        EXPECT_NEAR(std::get<3>(fresnel(1, 1.5)), 1. / 1.5, 1E-7);

        EXPECT_NEAR(std::get<0>(fresnel(-1, 1.5)), 0.04, 1E-7);
        EXPECT_NEAR(std::get<1>(fresnel(-1, 1.5)), 1, 1E-7);
        EXPECT_NEAR(std::get<2>(fresnel(-1, 1.5)), 1. / 1.5, 1E-7);
        EXPECT_NEAR(std::get<3>(fresnel(-1, 1.5)), 1.5, 1E-7);

        EXPECT_NEAR(std::get<0>(fresnel(1, 1. / 1.5)), 0.04, 1E-7);
        EXPECT_NEAR(std::get<1>(fresnel(1, 1. / 1.5)), -1, 1E-7);
        EXPECT_NEAR(std::get<2>(fresnel(1, 1. / 1.5)), 1. / 1.5, 1E-7);
        EXPECT_NEAR(std::get<3>(fresnel(1, 1. / 1.5)), 1.5, 1E-7);

        EXPECT_NEAR(std::get<0>(fresnel(-1, 1. / 1.5)), 0.04, 1E-7);
        EXPECT_NEAR(std::get<1>(fresnel(-1, 1. / 1.5)), 1, 1E-7);
        EXPECT_NEAR(std::get<2>(fresnel(-1, 1. / 1.5)), 1.5, 1E-7);
        EXPECT_NEAR(std::get<3>(fresnel(-1, 1. / 1.5)), 1. / 1.5, 1E-7);

        EXPECT_NEAR(std::get<0>(fresnel(0, 1.5)), 1, 1E-7);
        EXPECT_NEAR(std::get<1>(fresnel(0, 1.5)), ct_crit, 1E-7);
        EXPECT_NEAR(std::get<2>(fresnel(0, 1.5)), 1.5, 1E-7);
        EXPECT_NEAR(std::get<3>(fresnel(0, 1.5)), 1. / 1.5, 1E-7);

        EXPECT_NEAR(std::get<0>(fresnel(0, 1. / 1.5)), 1, 1E-7);
        EXPECT_NEAR(std::get<1>(fresnel(0, 1. / 1.5)), 0, 1E-7);
        EXPECT_NEAR(std::get<2>(fresnel(0, 1. / 1.5)), 1. / 1.5, 1E-7);
        EXPECT_NEAR(std::get<3>(fresnel(0, 1. / 1.5)), 1.5, 1E-7);

        // Spot check at 45 deg(1->1.5)
        // Source : http: // hyperphysics.phy-astr.gsu.edu/hbase/phyopt/freseq.html
        {
            auto [F, cos_theta_t, _nn, scale] = fresnel(std::cos(45 * EIGEN_PI / 180), 1.5);
            double cos_theta_t_ref            = -std::cos(28.1255057020557 * EIGEN_PI / 180);
            double F_ref                      = 0.5 * (0.09201336304552442 * 0.09201336304552442 + 0.3033370452904235 * 0.3033370452904235);
            double L                          = pow(scale * std::sqrt(1 - pow(std::cos(45 * EIGEN_PI / 180), 2)), 2) + cos_theta_t * cos_theta_t;
            EXPECT_NEAR(L, 1, 1E-7);
            EXPECT_NEAR(cos_theta_t, cos_theta_t_ref, 1E-7);
            EXPECT_NEAR(F, F_ref, 1E-7);
        }

        // 1.5->1
        {
            auto [F, cos_theta_t, _nn1, _nn2] = fresnel(std::cos(45 * EIGEN_PI / 180), 1 / 1.5);
            EXPECT_NEAR(F, 1, 1E-7);
            EXPECT_NEAR(cos_theta_t, 0, 1E-7);
        }

        {
            auto [F, cos_theta_t, _nn, scale] = fresnel(std::cos(10 * EIGEN_PI / 180), 1 / 1.5);
            double cos_theta_t_ref            = -std::cos(15.098086605159006 * EIGEN_PI / 180);
            double F_ref                      = 0.5 * (0.19046797197779405 * 0.19046797197779405 + 0.20949431963852014 * 0.20949431963852014);
            double L                          = pow(scale * std::sqrt(1 - pow(std::cos(10 * EIGEN_PI / 180), 2)), 2) + cos_theta_t * cos_theta_t;

            EXPECT_NEAR(L, 1, 1E-7);
            EXPECT_NEAR(cos_theta_t, cos_theta_t_ref, 1E-7);
            EXPECT_NEAR(F, F_ref, 1E-7);
        }
    }
}
