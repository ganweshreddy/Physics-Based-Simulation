#include "ocean.hpp"
#include "physsim_window.hpp"
#include "simulation.hpp"

#include <imgui.h>
#include <vislab/core/array.hpp>
#include <vislab/graphics/const_texture.hpp>
#include <vislab/graphics/diffuse_bsdf.hpp>
#include <vislab/graphics/mesh.hpp>
#include <vislab/graphics/perspective_camera.hpp>
#include <vislab/graphics/point_light.hpp>
#include <vislab/graphics/rectangle.hpp>
#include <vislab/graphics/scene.hpp>
#include <vislab/graphics/triangle.hpp>

#include <random>

using namespace vislab;

namespace physsim
{
    /**
     * @brief Ocean simulation using a simple Gerstner model.
     */
    class OceanSimulation : public Simulation
    {
    public:
        /**
         * @brief Initializes the scene.
         */
        void init() override
        {
            // grid resolution for the ocean surface
            Eigen::Vector2i resolution(100, 100);

            // allocate ocean simulation
            mOcean = std::make_shared<Ocean>(resolution, Eigen::Vector3d(-1, 1, 0.2), Eigen::Vector3d(2, 0, 0), Eigen::Vector3d(0, -2, 0));

            // add random waves
            std::default_random_engine rng;
            std::uniform_real_distribution<double> rnd;
            for (int i = 0; i < 50; ++i)
            {
                double lambda    = rnd(rng) * 1 + 0.1;
                double amplitude = rnd(rng) * 0.007;
                double angle     = rnd(rng) * 0.8 - 0.4;
                double steepness = rnd(rng) * 0.8;
                mOcean->waves.push_back(Wave{
                    lambda,    // Wavelength lambda in [m]
                    amplitude, // Amplitude in [m]
                    angle,     // direction angle in [rad]
                    steepness, // Steepness Q in [0,1]
                    0.0        // Phase offset
                });
            }

            // create triangle shape
            auto shape  = std::make_shared<Triangle>();
            shape->bsdf = std::make_shared<DiffuseBSDF>(std::make_shared<ConstTexture>(Spectrum(0.5, 1, 1)));
            shape->mesh = mOcean->mesh;
            scene->shapes.push_back(shape);

            // create ground plane
            auto ground  = std::make_shared<Rectangle>();
            ground->bsdf = std::make_shared<DiffuseBSDF>(std::make_shared<ConstTexture>(Spectrum(1, 1, 1)));
            scene->shapes.push_back(ground);

            // create a point light
            auto light       = std::make_shared<PointLight>();
            light->intensity = Eigen::Vector3d(100.0, 100.0, 100.0);
            light->transform.setMatrix(Eigen::Matrix4d::translate(Eigen::Vector3d(1, -5, 5)));

            // create a camera
            auto camera = std::make_shared<PerspectiveCamera>();
            camera->setLookAt(Eigen::Vector3d(0, 0, 0));
            camera->setPosition(Eigen::Vector3d(0, -2, 2));
            camera->setUp(Eigen::Vector3d(0, 0, 1));
            camera->setNear(0.01);
            camera->setFar(100);

            // add elements to scene
            scene->camera = camera;
            scene->lights.push_back(light);
        }

        /**
         * @brief Restarts the simulation.
         */
        void restart() override
        {
        }

        /**
         * @brief Advances the simulation one time step forward.
         * @param elapsedTime Elapsed time in milliseconds during the last frame.
         * @param totalTime Total time in milliseconds since the beginning of the first frame.
         * @param timeStep Time step of the simulation. Restarts when resetting the simulation.
         */
        void advance(double elapsedTime, double totalTime, int64_t timeStep) override
        {
            // update ocean mesh
            mOcean->advance(timeStep * 16.667);
        }

        /**
         * @brief Adds graphical user interface elements with imgui.
         */
        void gui() override
        {
        }

    private:
        /**
         * @brief Ocean simualtion.
         */
        std::shared_ptr<Ocean> mOcean;
    };
}

int main()
{
    physsim::PhyssimWindow window(
        800,       // width
        600,       // height
        "7_ocean", // title
        std::make_shared<physsim::OceanSimulation>(),
        false // fullscreen
    );

    return window.run();
}
