#include <vislab/core/array.hpp>
#include <vislab/field/regular_field.hpp>

#include "Eigen/Eigen"
#include "gtest/gtest.h"

namespace vislab
{
    TEST(field, regular)
    {
        // create a test grid
        auto grid = std::make_shared<RegularGrid3d>();
        grid->setDomain(Eigen::AlignedBox3d(Eigen::Vector3d(0.1, 0.2, 0.3), Eigen::Vector3d(1.1, 1.2, 1.3)));
        grid->setResolution(Eigen::Vector3i(16, 16, 16));

        // create a cubic scalar test data at grid coordinates
        auto array = std::make_shared<Array1d>();
        array->setSize(grid->getResolution().prod());
        for (int iz = 0; iz < grid->getResolution().z(); ++iz)
            for (int iy = 0; iy < grid->getResolution().y(); ++iy)
                for (int ix = 0; ix < grid->getResolution().x(); ++ix)
                {
                    Eigen::Index linearIndex = grid->getLinearIndex({ ix, iy, iz });
                    Eigen::Vector3d coord    = grid->getCoordAt(linearIndex);
                    double x = coord.x(), y = coord.y(), z = coord.z();
                    double s = 2 * x * x * x + x * x * y + 2 * x * y * y + 3 * y * y * y - 2 * y + x + x * x * z + 3 * z * z * y - 2 * z * z * z + z * 2;
                    array->setValue(linearIndex, s);
                }

        // place data and grid in field data structure
        auto field = std::make_shared<RegularSteadyScalarField3d>();
        field->setGrid(grid);
        field->setArray(array);
        field->setAccuracy(3); // third-order accurate reconstruction
        field->setBoundaryBehavior(0, EBoundaryBehavior::Clamp);
        field->setBoundaryBehavior(1, EBoundaryBehavior::Clamp);
        field->setBoundaryBehavior(2, EBoundaryBehavior::Clamp);

        for (int iz = 0; iz < grid->getResolution().z(); ++iz)
            for (int iy = 0; iy < grid->getResolution().y(); ++iy)
                for (int ix = 0; ix < grid->getResolution().x(); ++ix)
                {
                    Eigen::Index linearIndex = grid->getLinearIndex({ ix, iy, iz });
                    Eigen::Vector3d coord    = grid->getCoordAt(linearIndex);
                    double x = coord.x(), y = coord.y(), z = coord.z();
                    double s   = 2 * x * x * x + x * x * y + 2 * x * y * y + 3 * y * y * y - 2 * y + x + x * x * z + 3 * z * z * y - 2 * z * z * z + z * 2;
                    double sx  = 2 * x * z + 2 * y * y + 2 * x * y + 6 * x * x + 1;
                    double sy  = 3 * z * z + 9 * y * y + 4 * x * y + x * x - 2;
                    double sz  = -6 * z * z + 6 * y * z + x * x + 2;
                    double sxx = 2 * z + 2 * y + 12 * x;
                    double sxy = 4 * y + 2 * x;
                    double sxz = 2 * x;
                    double syy = 18 * y + 4 * x;
                    double syz = 6 * z;
                    double szz = 6 * y - 12 * z;

                    EXPECT_NEAR(s, field->sample(coord).x(), 1E-8);
                    EXPECT_NEAR(sx, field->sample_dx(coord).x(), 1E-8);
                    EXPECT_NEAR(sy, field->sample_dy(coord).x(), 1E-8);
                    EXPECT_NEAR(sz, field->sample_dz(coord).x(), 1E-8);
                    EXPECT_NEAR(sxx, field->sample_dxx(coord).x(), 1E-8);
                    EXPECT_NEAR(sxy, field->sample_dxy(coord).x(), 1E-8);
                    EXPECT_NEAR(sxz, field->sample_dxz(coord).x(), 1E-8);
                    EXPECT_NEAR(syy, field->sample_dyy(coord).x(), 1E-8);
                    EXPECT_NEAR(syz, field->sample_dyz(coord).x(), 1E-8);
                    EXPECT_NEAR(szz, field->sample_dzz(coord).x(), 1E-8);
                }
    }
}
