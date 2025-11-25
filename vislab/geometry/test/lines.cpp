#include "vislab/geometry/lines.hpp"
#include "vislab/geometry/attributes.hpp"

#include "vislab/core/array.hpp"

#include "Eigen/Eigen"
#include "gtest/gtest.h"

namespace vislab
{
    TEST(geometry, lines)
    {
        auto linesA  = std::make_unique<Lines2f>();
        auto linesA1 = linesA->createLine();
        linesA1->getVertices()->append(Eigen::Vector2f(1.0f, 2.0f));
        linesA1->getVertices()->append(Eigen::Vector2f(1.1f, 2.1f));
        linesA1->getVertices()->append(Eigen::Vector2f(1.2f, 2.2f));
        auto lineA1Float = linesA1->getAttributes()->create<Array1f>("float");
        lineA1Float->setSize(3);
        lineA1Float->setValue(0, 6.0f);
        lineA1Float->setValue(1, 6.1f);
        lineA1Float->setValue(2, 6.2f);
        EXPECT_NEAR(linesA1->length(), sqrtf(0.08f), 1E-6);

        auto linesA2 = linesA->createLine();
        linesA2->getVertices()->append(Eigen::Vector2f(3.0f, 4.0f));
        linesA2->getVertices()->append(Eigen::Vector2f(3.1f, 4.1f));
        linesA2->getVertices()->append(Eigen::Vector2f(3.2f, 4.2f));
        auto lineA2Float = linesA2->getAttributes()->create<Array1f>("float");
        lineA2Float->setSize(3);
        lineA2Float->setValue(0, 7.0f);
        lineA2Float->setValue(1, 7.1f);
        lineA2Float->setValue(2, 7.2f);
        EXPECT_NEAR(linesA2->length(), sqrtf(0.08f), 1E-6);

        auto linesB = linesA->clone();
        EXPECT_TRUE(linesB->isEqual(linesA.get()));
    }
}
