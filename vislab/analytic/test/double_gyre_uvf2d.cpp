#include "vislab/analytic/double_gyre_uvf2d.hpp"

#include "Eigen/Eigen"
#include "gtest/gtest.h"

namespace vislab
{
    void compare(const Eigen::Vector2d& a, const Eigen::Vector2d& b)
    {
        EXPECT_NEAR(a.x(), b.x(), 1E-10);
        EXPECT_NEAR(a.y(), b.y(), 1E-10);
    }

    TEST(analytic, double_gyre_uvf2d)
    {
        double A = 0.5, eps = 0.6, omega = 0.4;
        auto field = std::make_unique<DoubleGyreUVF2d>(A, eps, omega);
        Eigen::Vector3d pos(0.1, 0.2, 0.3);
        compare(field->sample(pos), Eigen::Vector2d(-0.3405365200287267, 0.7745178496335721));
        compare(field->sample_dx(pos), Eigen::Vector2d(-3.349039447371539, -0.4614968589102519));
        compare(field->sample_dy(pos), Eigen::Vector2d(0.7772748346149809, 3.349039447371539));
        compare(field->sample_dt(pos), Eigen::Vector2d(0.1741312479357944, -0.350871447472668));
        compare(field->sample_dxx(pos), Eigen::Vector2d(1.995526876571931, -6.087010272950341));
        compare(field->sample_dxy(pos), Eigen::Vector2d(7.644184777465769, -1.995526876571931));
        compare(field->sample_dyy(pos), Eigen::Vector2d(3.360960736807176, -7.644184777465769));
        compare(field->sample_dxt(pos), Eigen::Vector2d(1.517179131634284, 1.310815348478929));
        compare(field->sample_dyt(pos), Eigen::Vector2d(-0.3974546898205729, -1.517179131634284));
        compare(field->sample_dtt(pos), Eigen::Vector2d(-0.001510159260447787, -0.02892860253897779));
        compare(field->sample_dxxx(pos), Eigen::Vector2d(26.32042312557551, -1.375700143976786));
        compare(field->sample_dxxy(pos), Eigen::Vector2d(-4.554791429789536, -26.32042312557551));
        compare(field->sample_dxyy(pos), Eigen::Vector2d(33.0536944692, 4.554791429789536));
        compare(field->sample_dyyy(pos), Eigen::Vector2d(-7.671395128572018, -33.0536944692));
        compare(field->sample_dxxt(pos), Eigen::Vector2d(-5.66800663451839, 7.659976691265276));
        compare(field->sample_dxyt(pos), Eigen::Vector2d(-3.462962382192837, 5.66800663451839));
        compare(field->sample_dyyt(pos), Eigen::Vector2d(-1.718606530994299, 3.462962382192837));
        compare(field->sample_dxtt(pos), Eigen::Vector2d(0.1250881837083611, -0.9000130337641111));
        compare(field->sample_dytt(pos), Eigen::Vector2d(0.003446939521516861, -0.1250881837083611));
        compare(field->sample_dttt(pos), Eigen::Vector2d(-0.03238011990901046, 0.08530617854906686));
    }
}
