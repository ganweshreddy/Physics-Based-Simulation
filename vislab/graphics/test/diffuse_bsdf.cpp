#include "vislab/graphics/diffuse_bsdf.hpp"

#include "Eigen/Eigen"
#include "gtest/gtest.h"

namespace vislab
{
    TEST(graphics, diffuse_bsdf)
    {
        SurfaceInteraction si;
        si.position = Eigen::Vector3d(0, 0, 0);
        si.normal   = Eigen::Vector3d(0, 0, 1);
        si.wi       = Eigen::Vector3d(0, 0, 1);
        si.sh_frame = Frame(si.normal);

        DiffuseBSDF bsdf;

        for (int i = 0; i < 20; ++i)
        {
            double theta = i / 19.0 * (EIGEN_PI / 2);
            Eigen::Vector3d wo(std::sin(theta), 0, std::cos(theta));

            double v_pdf = bsdf.pdf(si, wo);
            double v_eval = bsdf.evaluate(si, wo).x();
                    
            EXPECT_NEAR(v_pdf, wo.z() / EIGEN_PI, 1E-7);
            EXPECT_NEAR(v_eval, 0.5 * wo.z() / EIGEN_PI, 1E-7);
        }
    }
}
