#include "vislab/geometry/surfaces.hpp"
#include "vislab/geometry/attributes.hpp"

#include "vislab/core/array.hpp"

#include "Eigen/Eigen"
#include "gtest/gtest.h"

namespace vislab
{
    TEST(geometry, surfaces)
    {
        auto surfacesA  = std::make_unique<Surfaces3f>();
        auto surfacesA1 = surfacesA->createSurface();
        surfacesA1->getVertices()->append(Eigen::Vector3f(1.0f, 2.0f, 3.2f));
        surfacesA1->getVertices()->append(Eigen::Vector3f(1.1f, 2.1f, 3.1f));
        surfacesA1->getVertices()->append(Eigen::Vector3f(1.2f, 2.2f, 3.0f));
        surfacesA1->recomputeBoundingBox();
        EXPECT_EQ(surfacesA1->getBoundingBox().min().x(), 1.0f);
        EXPECT_EQ(surfacesA1->getBoundingBox().min().y(), 2.0f);
        EXPECT_EQ(surfacesA1->getBoundingBox().min().z(), 3.0f);
        EXPECT_EQ(surfacesA1->getBoundingBox().max().x(), 1.2f);
        EXPECT_EQ(surfacesA1->getBoundingBox().max().y(), 2.2f);
        EXPECT_EQ(surfacesA1->getBoundingBox().max().z(), 3.2f);
        auto surfaceA1Float = surfacesA1->getAttributes()->create<Array1f>("float");
        surfaceA1Float->setSize(3);
        surfaceA1Float->setValue(0, 6.0f);
        surfaceA1Float->setValue(1, 6.1f);
        surfaceA1Float->setValue(2, 6.2f);

        auto surfacesA2 = surfacesA->createSurface();
        surfacesA2->getVertices()->append(Eigen::Vector3f(3.0f, 4.0f, 5.0f));
        surfacesA2->getVertices()->append(Eigen::Vector3f(3.1f, 4.1f, 5.2f));
        surfacesA2->getVertices()->append(Eigen::Vector3f(3.2f, 4.2f, 5.1f));
        surfacesA2->recomputeBoundingBox();
        EXPECT_EQ(surfacesA2->getBoundingBox().min().x(), 3.0f);
        EXPECT_EQ(surfacesA2->getBoundingBox().min().y(), 4.0f);
        EXPECT_EQ(surfacesA2->getBoundingBox().min().z(), 5.0f);
        EXPECT_EQ(surfacesA2->getBoundingBox().max().x(), 3.2f);
        EXPECT_EQ(surfacesA2->getBoundingBox().max().y(), 4.2f);
        EXPECT_EQ(surfacesA2->getBoundingBox().max().z(), 5.2f);
        auto surfaceA2Float = surfacesA2->getAttributes()->create<Array1f>("float");
        surfaceA2Float->setSize(3);
        surfaceA2Float->setValue(0, 7.0f);
        surfaceA2Float->setValue(1, 7.1f);
        surfaceA2Float->setValue(2, 7.2f);

        auto surfacesB = surfacesA->clone();
        EXPECT_TRUE(surfacesB->isEqual(surfacesA.get()));
        EXPECT_FALSE(surfaceA1Float->isEqual(surfaceA2Float.get()));
    }
}
