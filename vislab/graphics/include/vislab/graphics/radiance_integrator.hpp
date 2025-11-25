#pragma once

#include "image_fwd.hpp"

#include <vislab/core/object.hpp>

namespace vislab
{
    class Scene;

    /**
     * @brief Base class for radiance integrators
     */
    class RadianceIntegrator : public Object
    {
        VISLAB_INTERFACE(RadianceIntegrator, Object)

    public:
        /**
         * @brief Renders a scene onto a given sensor.
         * @param scene Scene to render.
         * @param camera Camera to render scene for.
         * @param image Image to render into.
         * @return True if successful.
         */
        virtual bool render(Scene* scene, Image3d* image) = 0;
    };
}
