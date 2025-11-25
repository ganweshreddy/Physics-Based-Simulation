#include "ocean.hpp"

#include <vislab/core/array.hpp>
#include <vislab/graphics/mesh.hpp>

namespace physsim
{
    Ocean::Ocean(const Eigen::Vector2i& resolution, const Eigen::Vector3d& origin, const Eigen::Vector3d& axis1, const Eigen::Vector3d& axis2)
        : mResolution(resolution)
        , mOrigin(origin)
        , mAxis1(axis1)
        , mAxis2(axis2)
    {
        // allocate mesh
        mesh = std::make_shared<vislab::Mesh>();

        // create buffer to hold the vertex positions
        mesh->positions = std::make_shared<vislab::Array3f>();
        mesh->positions->setSize(resolution.prod());

        // create index buffer
        mesh->indices = std::make_shared<vislab::Array3ui>();
        mesh->indices->setSize((resolution.x() - 1) * (resolution.y() - 1) * 2);
        int index = 0;
        for (int iy = 0; iy < resolution.y() - 1; ++iy)
            for (int ix = 0; ix < resolution.x() - 1; ++ix)
            {
                int i00 = (iy + 0) * resolution.x() + (ix + 0);
                int i01 = (iy + 1) * resolution.x() + (ix + 0);
                int i10 = (iy + 0) * resolution.x() + (ix + 1);
                int i11 = (iy + 1) * resolution.x() + (ix + 1);
                mesh->indices->setValue(index++, Eigen::Vector3u(i00, i01, i10));
                mesh->indices->setValue(index++, Eigen::Vector3u(i01, i11, i10));
            }
    }

    void Ocean::advance(double totalTime)
    {
#ifndef _DEBUG
#pragma omp parallel for
#endif
        for (int64_t i = 0; i < mResolution.prod(); ++i)
        {
            // vertex grid indices
            int ix = i % mResolution.x();
            int iy = i / mResolution.x();

            double x = ix / (mResolution.x() - 1.);
            double y = iy / (mResolution.y() - 1.);

            Eigen::Vector3d pos = mOrigin + x * mAxis1 + y * mAxis2;

            Eigen::Vector3d surfacePosition = position(pos, totalTime / 1000); // conversion from milliseconds to seconds
            mesh->positions->setValue(i, surfacePosition.cast<float>());
        }

        // recompute the vertex normals
        mesh->recomputeVertexNormals();

        // raise event that geometry changed
        mesh->positionsChanged.notify(mesh.get());
    }

    const Eigen::Vector2i& Ocean::resolution() const
    {
        return mResolution;
    }

    Eigen::Vector3d Ocean::position(const Eigen::Vector3d& pos, double t) const
    {
        Eigen::Vector3d result = pos;
        for (auto& wave : waves)
            result += wave.offset(pos.xy(), t);
        return result;
    }
}
