#include "vislab/graphics/ray_box_intersect.hpp"
#include "vislab/graphics/ray.hpp"

#include "Eigen/Eigen"
#include "gtest/gtest.h"

namespace vislab
{
    TEST(graphics, ray_box_intersect)
    {
        double EPS = 1E-10;

        // generate a test box
        Eigen::AlignedBox3d box(Eigen::Vector3d(-1, -2, -3), Eigen::Vector3d(2, 3, 4));

        // from inside (one hit)
        Ray3d ray1(Eigen::Vector3d(0, 1, 2), Eigen::Vector3d(1, 2, 3));
        Ray3d ray1n(Eigen::Vector3d(0, 1, 2), Eigen::Vector3d(1, 2, 3).normalized());

        // from outside (two hits)
        Ray3d ray2(Eigen::Vector3d(4, 5, 6), Eigen::Vector3d(-4, -5, -6));
        Ray3d ray2n(Eigen::Vector3d(4, 5, 6), Eigen::Vector3d(-4, -5, -6).normalized());

        // shoot at corner
        Ray3d ray3(Eigen::Vector3d(-2, -3, -4), Eigen::Vector3d(1, 1, 1));
        Ray3d ray3n(Eigen::Vector3d(-2, -3, -4), Eigen::Vector3d(1, 1, 1).normalized());

        // axis-parallel hit
        Ray3d ray4(Eigen::Vector3d(-3, 0, 0), Eigen::Vector3d(2, 0, 0));
        Ray3d ray4n(Eigen::Vector3d(-3, 0, 0), Eigen::Vector3d(2, 0, 0).normalized());

        // miss
        Ray3d ray5(Eigen::Vector3d(5, 6, 7), Eigen::Vector3d(1, 2, 3));
        Ray3d ray5n(Eigen::Vector3d(5, 6, 7), Eigen::Vector3d(1, 2, 3).normalized());

        // any hits with non-normalized ray
        EXPECT_TRUE(RayBoxIntersection::anyHit(ray1, box));
        EXPECT_TRUE(RayBoxIntersection::anyHit(ray2, box));
        EXPECT_TRUE(RayBoxIntersection::anyHit(ray3, box));
        EXPECT_TRUE(RayBoxIntersection::anyHit(ray4, box));
        EXPECT_FALSE(RayBoxIntersection::anyHit(ray5, box));

        // any hits with normalized ray
        EXPECT_TRUE(RayBoxIntersection::anyHit(ray1n, box));
        EXPECT_TRUE(RayBoxIntersection::anyHit(ray2n, box));
        EXPECT_TRUE(RayBoxIntersection::anyHit(ray3n, box));
        EXPECT_TRUE(RayBoxIntersection::anyHit(ray4n, box));
        EXPECT_FALSE(RayBoxIntersection::anyHit(ray5n, box));

        // closest hits with non-normalized ray
        Eigen::Vector2d t1, t2, t3, t4, t5;
        EXPECT_TRUE(RayBoxIntersection::closestHit(ray1, box, t1));
        EXPECT_NEAR(t1.x(), -1, EPS);
        EXPECT_NEAR(t1.y(), 0.66666666666666666, EPS);

        EXPECT_TRUE(RayBoxIntersection::closestHit(ray2, box, t2));
        EXPECT_NEAR(t2.x(), 0.5, EPS);
        EXPECT_NEAR(t2.y(), 1.25, EPS);

        EXPECT_TRUE(RayBoxIntersection::closestHit(ray3, box, t3));
        EXPECT_NEAR(t3.x(), 1, EPS);
        EXPECT_NEAR(t3.y(), 4, EPS);

        EXPECT_TRUE(RayBoxIntersection::closestHit(ray4, box, t4));
        EXPECT_NEAR(t4.x(), 1, EPS);
        EXPECT_NEAR(t4.y(), 2.5, EPS);

        EXPECT_FALSE(RayBoxIntersection::closestHit(ray5, box, t5));
        EXPECT_NEAR(t5.x(), -3.3333333333333333, EPS);
        EXPECT_NEAR(t5.y(), -3, EPS);

        // closest hits with normalized ray
        Eigen::Vector2d t1n, t2n, t3n, t4n, t5n;
        EXPECT_TRUE(RayBoxIntersection::closestHit(ray1n, box, t1n));
        EXPECT_NEAR(t1n.x(), -3.7416573867739413, EPS);
        EXPECT_NEAR(t1n.y(), 2.4944382578492941, EPS);

        EXPECT_TRUE(RayBoxIntersection::closestHit(ray2n, box, t2n));
        EXPECT_NEAR(t2n.x(), 4.3874821936960613, EPS);
        EXPECT_NEAR(t2n.y(), 10.968705484240152, EPS);

        EXPECT_TRUE(RayBoxIntersection::closestHit(ray3n, box, t3n));
        EXPECT_NEAR(t3n.x(), 1.7320508075688770, EPS);
        EXPECT_NEAR(t3n.y(), 6.9282032302755079, EPS);

        EXPECT_TRUE(RayBoxIntersection::closestHit(ray4n, box, t4n));
        EXPECT_NEAR(t4n.x(), 2, EPS);
        EXPECT_NEAR(t4n.y(), 5, EPS);

        EXPECT_FALSE(RayBoxIntersection::closestHit(ray5n, box, t5n));
        EXPECT_NEAR(t5n.x(), -12.472191289246471, EPS);
        EXPECT_NEAR(t5n.y(), -11.224972160321824, EPS);
    }
}
