#pragma once

#include "interaction.hpp"

#include "frame.hpp"
#include "math.hpp"
#include "ray.hpp"

namespace vislab
{
    class PositionSample;
    class Scene;
    class Light;
    class Shape;
    class BSDF;
    class RayDifferential3d;

    /**
     * @brief Represents the interaction on a surface.
     */
    class SurfaceInteraction : public Interaction
    {
    public:
        /**
         * @brief Empty constructor.
         */
        SurfaceInteraction();

        /**
         * @brief Constructor.
         * @param ps
         */
        SurfaceInteraction(const PositionSample& ps);

        /**
         * @brief Incident ray direction in the local shading frame.
         */
        Eigen::Vector3d wi;

        /**
         * @brief Geometric normal at the interaction.
         */
        Eigen::Vector3d normal;

        /**
         * @brief Shape that was interacted with.
         */
        const Shape* shape;

        /**
         * @brief Index of the intersected primitive.
         */
        int64_t prim_index;

        /**
         * @brief UV coordinates of the surface.
         */
        Eigen::Vector2d uv;

        /**
         * @brief Shading frame.
         */
        Frame sh_frame;

        /**
         * @brief Position partials wrt. to u-parameterization.
         */
        Eigen::Vector3d dp_du;

        /**
         * @brief Position partials wrt. to v-parameterization.
         */
        Eigen::Vector3d dp_dv;

        /**
         * @brief UV partials wrt. changes in screen-space.
         */
        Eigen::Vector2d duv_dx;

        /**
         * @brief UV partials wrt. changes in screen-space.
         */
        Eigen::Vector2d duv_dy;

        /**
         * @brief Convert a local shading-space vector into world space.
         * @param v Local shading-space vector to transform.
         * @return Transformed vector in world space.
         */
        Eigen::Vector3d toWorld(const Eigen::Vector3d& v) const;

        /**
         * @brief Convert a world-space vector into local shading coordinates.
         * @param v World-space vector to transform.
         * @return Transformed vector in local shading coordinates.
         */
        Eigen::Vector3d toLocal(const Eigen::Vector3d& v) const;

        /**
         * @brief Returns the light associated with the intersection (if any).
         * @param scene Scene to fall back onto when getting an environment light.
         * @return Light associated with the intersection.
         */
        const Light* light(const Scene* scene) const;

        /**
         * @brief Initialize local shading frame using Gram-schmidt orthogonalization.
         */
        void initializeShadingFrame();

        /**
         * @brief Returns the BSDF of the intersected shape.
         * @param ray Differential ray that is used to perform texture filtering.
         * @return BSDF of intersected shape.
         */
        const BSDF* bsdf(const RayDifferential3d& ray);

        /**
         * @brief Returns the BSDF of the intersected shape.
         * @return BSDF of intersected shape.
         */
        const BSDF* bsdf() const;

        /**
         * @brief Computes UV partials with respect to screen space from a ray differential.
         * @param ray Ray differential to compute UV partials with respect to screen space from.
         */
        void computeUVPartials(const RayDifferential3d& ray);

        /**
         * @brief Checks is surface interaction has UV partials computed.
         */
        bool hasUVPartials() const;

        /**
         * @brief Spawn a semi-infinite ray towards the given direction.
         * @param direction Direction to trave.
         * @return Spawned ray.
         */
        Ray3d spawnRay(const Eigen::Vector3d& direction) const;

    private:
        /**
         * Compute an offset position, used when spawning a ray from this
         * interaction. When the interaction is on the surface of a shape, the
         * position is offset along the surface normal to prevent self intersection.
         */
        Eigen::Vector3d offsetPosition(const Eigen::Vector3d& direction) const;
    };
}
