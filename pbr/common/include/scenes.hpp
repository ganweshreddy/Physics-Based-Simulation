#pragma once

#include <vislab/graphics/scene.hpp>

namespace pbr
{
    /**
     * @brief Collection of common scenes.
     */
    class Scenes
    {
    public:
        /**
         * @brief Creates a scene that contains a cornell box, an area light, and a camera.
         * @param resolution Resolution of the viewport, which is used to compute the aspect ratio of the camera.
         * @return Cornell box scene.
         */
        static std::shared_ptr<vislab::Scene> create_box_scene(const Eigen::Vector2i& resolution, bool addLight = true);
    };
}
