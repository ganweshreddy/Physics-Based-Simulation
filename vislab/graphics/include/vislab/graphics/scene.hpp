#pragma once

#include "bounding_volume_hierarchy.hpp"
#include "direction_sample.hpp"
#include "ray_fwd.hpp"
#include "spectrum.hpp"
#include "surface_interaction.hpp"

#include <vislab/core/data.hpp>

#include <memory>
#include <vector>

namespace vislab
{
    class Shape;
    class Light;
    class Camera;

    /**
     * @brief Representation of the scene to render.
     */
    class Scene : public Data
    {
        VISLAB_OBJECT(Scene, Data)

    public:
        /**
         * @brief Computes the ray intersection with the shapes in the scene.
         * @param ray Ray to compute intersection with.
         * @return Surface interaction result.
         */
        SurfaceInteraction intersect(const Ray3d& ray) const;

        /**
         * @brief Determines if there is any hit with the shapes in the scene.
         * @param ray Ray to compute intersection with.
         * @return True if any object was hit.
         */
        bool anyHit(const Ray3d& ray) const;

        /**
         * @brief Given an arbitrary reference point in the scene, this method samples a direction from the reference point to towards an emitter.
         * @param it Interaction to start the ray towards the light from.
         * @param sample Uniformly distributed random number in [0,1]^2.
         * @param testVisibility True if a visibility test is performed.
         * @return Direction sample and evaluation from the light source.
         */
        std::pair<DirectionSample, Spectrum> sampleLightDirection(const Interaction& it, const Eigen::Vector2d& sample, bool testVisibility = true) const;

        /**
         * @brief Evaluate the probability density of the sampleLightDirection() technique given an filled-in DirectionSample record.
         * @param it Interaction to start the ray towards the light from.
         * @param ds Direction sample that was computed by sampleLightDirection.
         * @return Probability density for choosing the direction was sampled.
         */
        double pdfLightDirection(const Interaction& it, const DirectionSample& ds) const;

        /**
         * @brief Collection of shapes in the scene.
         */
        std::vector<std::shared_ptr<Shape>> shapes;

        /**
         * @brief Collection of lights in the scene.
         */
        std::vector<std::shared_ptr<Light>> lights;

        /**
         * @brief Environment light.
         */
        std::shared_ptr<Light> environment;

        /**
         * @brief Camera to render the scene from.
         */
        std::shared_ptr<Camera> camera;

        /**
         * @brief Rebuilds the ray traversal acceleration data structure.
         */
        void buildAccelerationTree();

    private:
        /**
         * @brief Bounding volume hierarchy for fast ray traversal.
         */
        std::shared_ptr<BoundingVolumeHierarchy<Shape, 3>> mBoundingVolumeHierarchy;
    };
}
