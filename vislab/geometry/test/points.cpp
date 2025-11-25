#include "vislab/geometry/points.hpp"
#include "vislab/geometry/attributes.hpp"

#include "vislab/core/array.hpp"

#include "Eigen/Eigen"
#include "gtest/gtest.h"

namespace vislab
{
    TEST(geometry, points)
    {
        auto pointsA = std::make_unique<Points2f>();
        pointsA->getVertices()->append(Eigen::Vector2f(1.0f, 2.0f));
        pointsA->getVertices()->append(Eigen::Vector2f(1.1f, 2.1f));
        pointsA->getVertices()->append(Eigen::Vector2f(1.2f, 2.2f));

        auto pointsArev = std::make_unique<Points2f>();
        pointsArev->getVertices()->setValues({ Eigen::Vector2f(1.2f, 2.2f), Eigen::Vector2f(1.1f, 2.1f), Eigen::Vector2f(1.0f, 2.0f) });

        auto pointsB = std::make_unique<Points2f>();
        pointsB->getVertices()->setValues({ Eigen::Vector2f(11.0f, 12.0f), Eigen::Vector2f(11.1f, 12.1f), Eigen::Vector2f(11.2f, 12.2f), Eigen::Vector2f(11.3f, 12.3f) });

        auto pointsAB_gt = std::make_unique<Points2f>();
        pointsAB_gt->getVertices()->setValues({ Eigen::Vector2f(1.0f, 2.0f), Eigen::Vector2f(1.1f, 2.1f), Eigen::Vector2f(1.2f, 2.2f), Eigen::Vector2f(11.0f, 12.0f), Eigen::Vector2f(11.1f, 12.1f), Eigen::Vector2f(11.2f, 12.2f), Eigen::Vector2f(11.3f, 12.3f) });

        auto pointsBA_gt = std::make_unique<Points2f>();
        pointsBA_gt->getVertices()->setValues({ Eigen::Vector2f(11.0f, 12.0f), Eigen::Vector2f(11.1f, 12.1f), Eigen::Vector2f(11.2f, 12.2f), Eigen::Vector2f(11.3f, 12.3f), Eigen::Vector2f(1.0f, 2.0f), Eigen::Vector2f(1.1f, 2.1f), Eigen::Vector2f(1.2f, 2.2f) });

        auto pointsTest = std::make_unique<Points2f>();
        pointsTest->getVertices()->setValues({ Eigen::Vector2f(1.0f, 2.0f), Eigen::Vector2f(1.1f, 2.1f), Eigen::Vector2f(1.2f, 2.2f) });

        EXPECT_TRUE(pointsA->isEqual(pointsTest.get()));
        EXPECT_FALSE(pointsA->isEqual(pointsB.get()));

        pointsTest->append(pointsB.get());
        EXPECT_TRUE(pointsTest->isEqual(pointsAB_gt.get()));

        pointsTest->removeLast(pointsB->getNumVertices());
        EXPECT_TRUE(pointsTest->isEqual(pointsA.get()));

        pointsTest->getVertices()->setValues({ Eigen::Vector2f(1.0f, 2.0f), Eigen::Vector2f(1.1f, 2.1f), Eigen::Vector2f(1.2f, 2.2f) });
        pointsTest->prepend(pointsB.get());
        EXPECT_TRUE(pointsTest->isEqual(pointsBA_gt.get()));

        pointsTest->removeFirst(pointsB->getNumVertices());
        EXPECT_TRUE(pointsTest->isEqual(pointsA.get()));

        pointsTest->getVertices()->setValues({ Eigen::Vector2f(1.0f, 2.0f), Eigen::Vector2f(1.1f, 2.1f), Eigen::Vector2f(1.2f, 2.2f) });
        pointsTest->reverse();
        EXPECT_TRUE(pointsTest->isEqual(pointsArev.get()));

        auto pointsAattr = std::make_unique<Points2f>();
        pointsAattr->getVertices()->setValues({ Eigen::Vector2f(1.0f, 2.0f), Eigen::Vector2f(1.1f, 2.1f), Eigen::Vector2f(1.2f, 2.2f) });
        auto pointsAattr_float = pointsAattr->getAttributes()->create<Array1f>("float");
        pointsAattr_float->setValues({ 30.f, 31.f, 32.f });
        auto pointsAattr_vec2f = pointsAattr->getAttributes()->create<Array2f>("Vec2f");
        pointsAattr_vec2f->setValues({ Eigen::Vector2f(40.f, 50.f), Eigen::Vector2f(41.f, 51.f), Eigen::Vector2f(42.f, 52.f) });

        auto pointsAattrrev = std::make_unique<Points2f>();
        pointsAattrrev->getVertices()->setValues({ Eigen::Vector2f(1.2f, 2.2f), Eigen::Vector2f(1.1f, 2.1f), Eigen::Vector2f(1.0f, 2.0f) });
        auto pointsAattrrev_float = pointsAattrrev->getAttributes()->create<Array1f>("float");
        pointsAattrrev_float->setValues({ 32.f, 31.f, 30.f });
        auto pointsAattrrev_vec2f = pointsAattrrev->getAttributes()->create<Array2f>("Vec2f");
        pointsAattrrev_vec2f->setValues({ Eigen::Vector2f(42.f, 52.f), Eigen::Vector2f(41.f, 51.f), Eigen::Vector2f(40.f, 50.f) });

        auto pointsBattr = std::make_unique<Points2f>();
        pointsBattr->getVertices()->setValues({ Eigen::Vector2f(11.0f, 12.0f), Eigen::Vector2f(11.1f, 12.1f), Eigen::Vector2f(11.2f, 12.2f), Eigen::Vector2f(11.3f, 12.3f) });
        auto pointsBattr_float = pointsBattr->getAttributes()->create<Array1f>("float");
        pointsBattr_float->setValues({ 130.f, 131.f, 132.f, 133.f });
        auto pointsBattr_vec2f = pointsBattr->getAttributes()->create<Array2f>("Vec2f");
        pointsBattr_vec2f->setValues({ Eigen::Vector2f(140.f, 150.f), Eigen::Vector2f(141.f, 151.f), Eigen::Vector2f(142.f, 152.f), Eigen::Vector2f(143.f, 153.f) });

        auto pointsAB_gtattr = std::make_unique<Points2f>();
        pointsAB_gtattr->getVertices()->setValues({ Eigen::Vector2f(1.0f, 2.0f), Eigen::Vector2f(1.1f, 2.1f), Eigen::Vector2f(1.2f, 2.2f), Eigen::Vector2f(11.0f, 12.0f), Eigen::Vector2f(11.1f, 12.1f), Eigen::Vector2f(11.2f, 12.2f), Eigen::Vector2f(11.3f, 12.3f) });
        auto pointsAB_gtattr_float = pointsAB_gtattr->getAttributes()->create<Array1f>("float");
        pointsAB_gtattr_float->setValues({ 30.f, 31.f, 32.f, 130.f, 131.f, 132.f, 133.f });
        auto pointsAB_gtattr_vec2f = pointsAB_gtattr->getAttributes()->create<Array2f>("Vec2f");
        pointsAB_gtattr_vec2f->setValues({ Eigen::Vector2f(40.f, 50.f), Eigen::Vector2f(41.f, 51.f), Eigen::Vector2f(42.f, 52.f), Eigen::Vector2f(140.f, 150.f), Eigen::Vector2f(141.f, 151.f), Eigen::Vector2f(142.f, 152.f), Eigen::Vector2f(143.f, 153.f) });

        auto pointsBA_gtattr = std::make_unique<Points2f>();
        pointsBA_gtattr->getVertices()->setValues({ Eigen::Vector2f(11.0f, 12.0f), Eigen::Vector2f(11.1f, 12.1f), Eigen::Vector2f(11.2f, 12.2f), Eigen::Vector2f(11.3f, 12.3f), Eigen::Vector2f(1.0f, 2.0f), Eigen::Vector2f(1.1f, 2.1f), Eigen::Vector2f(1.2f, 2.2f) });
        auto pointsBA_gtattr_float = pointsBA_gtattr->getAttributes()->create<Array1f>("float");
        pointsBA_gtattr_float->setValues({ 130.f, 131.f, 132.f, 133.f, 30.f, 31.f, 32.f });
        auto pointsBA_gtattr_vec2f = pointsBA_gtattr->getAttributes()->create<Array2f>("Vec2f");
        pointsBA_gtattr_vec2f->setValues({ Eigen::Vector2f(140.f, 150.f), Eigen::Vector2f(141.f, 151.f), Eigen::Vector2f(142.f, 152.f), Eigen::Vector2f(143.f, 153.f), Eigen::Vector2f(40.f, 50.f), Eigen::Vector2f(41.f, 51.f), Eigen::Vector2f(42.f, 52.f) });

        pointsTest = std::make_unique<Points2f>();
        pointsTest->getVertices()->setValues({ Eigen::Vector2f(1.0f, 2.0f), Eigen::Vector2f(1.1f, 2.1f), Eigen::Vector2f(1.2f, 2.2f) });
        auto pointsAattr_float_ = pointsTest->getAttributes()->create<Array1f>("float");
        pointsAattr_float_->setValues({ 30.f, 31.f, 32.f });
        auto pointsAattr_vec2f_ = pointsTest->getAttributes()->create<Array2f>("Vec2f");
        pointsAattr_vec2f_->setValues({ Eigen::Vector2f(40.f, 50.f), Eigen::Vector2f(41.f, 51.f), Eigen::Vector2f(42.f, 52.f) });
        pointsTest->append(pointsBattr.get());
        EXPECT_TRUE(pointsTest->isEqual(pointsAB_gtattr.get()));
        pointsTest->removeLast(pointsBattr->getNumVertices());
        EXPECT_TRUE(pointsTest->isEqual(pointsAattr.get()));

        pointsTest = std::make_unique<Points2f>();
        pointsTest->getVertices()->setValues({ Eigen::Vector2f(1.0f, 2.0f), Eigen::Vector2f(1.1f, 2.1f), Eigen::Vector2f(1.2f, 2.2f) });
        auto pointsAattr_float__ = pointsTest->getAttributes()->create<Array1f>("float");
        pointsAattr_float__->setValues({ 30.f, 31.f, 32.f });
        auto pointsAattr_vec2f__ = pointsTest->getAttributes()->create<Array2f>("Vec2f");
        pointsAattr_vec2f__->setValues({ Eigen::Vector2f(40.f, 50.f), Eigen::Vector2f(41.f, 51.f), Eigen::Vector2f(42.f, 52.f) });
        pointsTest->prepend(pointsBattr.get());
        EXPECT_TRUE(pointsTest->isEqual(pointsBA_gtattr.get()));
        pointsTest->removeFirst(pointsBattr->getNumVertices());
        EXPECT_TRUE(pointsTest->isEqual(pointsAattr.get()));

        pointsTest->reverse();
        EXPECT_TRUE(pointsTest->isEqual(pointsAattrrev.get()));
        pointsTest->reverse();

        pointsTest->getAttributes()->clear();
        EXPECT_TRUE(pointsTest->isEqual(pointsA.get()));

        pointsTest->setAttributes(std::shared_ptr<Attributes>(pointsAattr->getAttributes()->clone()));
        EXPECT_TRUE(pointsTest->isEqual(pointsAattr.get()));

        auto cloned = pointsTest->clone();
        EXPECT_TRUE(cloned->isEqual(pointsTest.get()));
    }
}
