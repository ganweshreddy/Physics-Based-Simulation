#include "vislab/graphics/warp.hpp"

#include "Eigen/Eigen"
#include "gtest/gtest.h"

namespace vislab
{
    TEST(graphics, warp)
    {
        EXPECT_NEAR((Warp::square_to_cosine_hemisphere(Eigen::Vector2d(0.5, 0.5)) - Eigen::Vector3d(0, 0, 1)).squaredNorm(), 0, 1E-7);
        EXPECT_NEAR((Warp::square_to_cosine_hemisphere(Eigen::Vector2d(0.5, 0.0)) - Eigen::Vector3d(0, -1, 0)).squaredNorm(), 0, 1E-7);

        EXPECT_NEAR((Warp::square_to_uniform_cone(Eigen::Vector2d(0.5, 0.5), 1) - Eigen::Vector3d(0, 0, 1)).squaredNorm(), 0, 1E-7);
        EXPECT_NEAR((Warp::square_to_uniform_cone(Eigen::Vector2d(0.5, 0.0), 1) - Eigen::Vector3d(0, 0, 1)).squaredNorm(), 0, 1E-7);
        EXPECT_NEAR((Warp::square_to_uniform_cone(Eigen::Vector2d(0.5, 0.0), 0) - Eigen::Vector3d(0, -1, 0)).squaredNorm(), 0, 1E-7);

        EXPECT_NEAR((Warp::square_to_uniform_disk_concentric(Eigen::Vector2d(0.0, 0.0)) - Eigen::Vector2d(-1. / sqrt(2), -1. / sqrt(2))).squaredNorm(), 0, 1E-7);
        EXPECT_NEAR((Warp::square_to_uniform_disk_concentric(Eigen::Vector2d(0.5, 0.5)) - Eigen::Vector2d(0, 0)).squaredNorm(), 0, 1E-7);

        EXPECT_NEAR((Warp::square_to_uniform_sphere(Eigen::Vector2d(0, 0)) - Eigen::Vector3d(0, 0, 1)).squaredNorm(), 0, 1E-7);
        EXPECT_NEAR((Warp::square_to_uniform_sphere(Eigen::Vector2d(0, 1)) - Eigen::Vector3d(0, 0, -1)).squaredNorm(), 0, 1E-7);
        EXPECT_NEAR((Warp::square_to_uniform_sphere(Eigen::Vector2d(0.5, 0.5)) - Eigen::Vector3d(-1, 0, 0)).squaredNorm(), 0, 1E-7);
    }
}
