#include "physsim_window.hpp"
#include "simulation.hpp"

#include <imgui.h>
#include <vislab/graphics/const_texture.hpp>
#include <vislab/graphics/diffuse_bsdf.hpp>
#include <vislab/graphics/perspective_camera.hpp>
#include <vislab/graphics/point_light.hpp>
#include <vislab/graphics/scene.hpp>
#include <vislab/graphics/sphere.hpp>

using namespace vislab;

namespace physsim
{
    /**
     * @brief Basic example for a physsim simulation.
     */
    class EmptySimulation : public Simulation
    {
    public:
        /**
         * @brief Initializes the scene.
         */
        void init() override
        {
            // create texture with constant color
            mTex        = std::make_shared<ConstTexture>();
            mTex->color = Spectrum(0.0, 1.0, 1.0);

            // create a sphere and assign the bsdf
            auto sphere  = std::make_shared<Sphere>();
            sphere->bsdf = std::make_shared<DiffuseBSDF>(mTex);

            // create a point light
            auto light       = std::make_shared<PointLight>();
            light->intensity = Eigen::Vector3d(50.0, 50.0, 50.0);
            light->transform.setMatrix(Eigen::Matrix4d::translate(Eigen::Vector3d(-3, -3, 3)));

            // create a camera
            auto camera = std::make_shared<PerspectiveCamera>();
            camera->setLookAt(Eigen::Vector3d(0, 0, 0));
            camera->setPosition(Eigen::Vector3d(-3.5, 0, 0));
            camera->setUp(Eigen::Vector3d(0, 0, 1));
            camera->setNear(0.01);
            camera->setFar(10);

            // add elements to scene
            scene->camera = camera;
            scene->shapes.push_back(sphere);
            scene->lights.push_back(light);
        }

        /**
         * @brief Restarts the simulation.
         */
        void restart() override
        {
            mTex->color = Spectrum(0, 1, 1);
        }

        /**
         * @brief Advances the simulation one time step forward.
         * @param elapsedTime Elapsed time in milliseconds during the last frame.
         * @param totalTime Total time in milliseconds since the beginning of the first frame.
         * @param timeStep Time step of the simulation. Restarts when resetting the simulation.
         */
        void advance(double elapsedTime, double totalTime, int64_t timeStep) override
        {
            // compute time variable since start of simulation from the timestep.
            const double time = timeStep * 0.02;

            // animate color of the texture over time
            mTex->color = Spectrum(
                std::sin(time),
                std::cos(time),
                1);
        }

        /**
         * @brief Adds graphical user interface elements with imgui.
         */
        void gui() override
        {
        }

    private:
        /**
         * @brief Constant texture used to define the color of the sphere.
         */
        std::shared_ptr<ConstTexture> mTex;
    };
}

int main()
{
    physsim::PhyssimWindow window(
        800,       // width
        600,       // height
        "0_empty", // title
        std::make_shared<physsim::EmptySimulation>(),
        false // fullscreen
    );

    return window.run();
}
