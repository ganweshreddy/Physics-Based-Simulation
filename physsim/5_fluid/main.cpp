#include "fluid_solver.hpp"
#include "physsim_window.hpp"
#include "simulation.hpp"

#include <imgui.h>
#include <vislab/graphics/colormap_texture.hpp>
#include <vislab/graphics/diffuse_bsdf.hpp>
#include <vislab/graphics/perspective_camera.hpp>
#include <vislab/graphics/rectangle.hpp>
#include <vislab/graphics/scene.hpp>

using namespace vislab;

namespace physsim
{
    /**
     * @brief Fluid simulation.
     */
    class FluidSimulation : public Simulation
    {
    public:
        /**
         * @brief Initializes the scene.
         */
        void init() override
        {
            // simulation settings
            mBoundary       = FluidSolver::EBoundary::Closed;
            mPressureSolver = FluidSolver::EPressureSolver::Iterative;
            mAccuracy       = 1E-3;
            mIterations     = 1000;

            // allocate a fluid solver
            Eigen::AlignedBox2d domain(
                Eigen::Vector2d(0.0, 0.0),
                Eigen::Vector2d(1.0, 1.5));
            Eigen::Vector2i resolution(128, 192);
            mFluidSolver               = std::make_shared<FluidSolver>(domain, resolution);
            mFluidSolver->sourceRegion = Eigen::AlignedBox2d(
                Eigen::Vector2d(0.45, 0.1),
                Eigen::Vector2d(0.55, 0.15));

            // create plane to show the field
            auto plane  = std::make_shared<Rectangle>();
            plane->bsdf = std::make_shared<DiffuseBSDF>(mFluidSolver->densityTexture);
            plane->transform.setMatrix(Eigen::Matrix4d::scale(Eigen::Vector3d(1., 1.5, 1.)));
            scene->shapes.push_back(plane);

            // create a camera
            auto camera = std::make_shared<PerspectiveCamera>();
            camera->setLookAt(Eigen::Vector3d(0, 0, 0));
            camera->setPosition(Eigen::Vector3d(0, 0, -4));
            camera->setUp(Eigen::Vector3d(0, 1, 0));
            camera->setNear(0.01);
            camera->setFar(100);

            // add elements to scene
            scene->camera = camera;
        }

        /**
         * @brief Restarts the simulation.
         */
        void restart() override
        {
            mFluidSolver->reset();
        }

        /**
         * @brief Advances the simulation one time step forward.
         * @param elapsedTime Elapsed time in milliseconds during the last frame.
         * @param totalTime Total time in milliseconds since the beginning of the first frame.
         * @param timeStep Time step of the simulation. Restarts when resetting the simulation.
         */
        void advance(double elapsedTime, double totalTime, int64_t timeStep) override
        {
            mFluidSolver->iterations     = mIterations;
            mFluidSolver->accuracy       = mAccuracy;
            mFluidSolver->boundary       = mBoundary;
            mFluidSolver->pressureSolver = mPressureSolver;
            mFluidSolver->advance();
        }

        /**
         * @brief Adds graphical user interface elements with imgui.
         */
        void gui() override
        {
            ImGui::PushItemWidth(100);

            ImGui::SliderInt("iterations", &mIterations, 1, 1000);
            ImGui::SliderFloat("accuracy", &mAccuracy, 1E-4, 1E-2);
            ImGui::Combo("boundary", (int*)&mBoundary, "open\0closed\0\0");
            ImGui::Combo("solver", (int*)&mPressureSolver, "iterative\0linear\0\0");

            ImGui::PopItemWidth();
        }

    private:
        /**
         * @brief Fluid solver.
         */
        std::shared_ptr<FluidSolver> mFluidSolver;

        /**
         * @brief Boundary setting.
         */
        FluidSolver::EBoundary mBoundary;

        /**
         * @brief Pressure solver to use.
         */
        FluidSolver::EPressureSolver mPressureSolver;

        /**
         * @brief Accuracy of the iterative pressure solver.
         */
        float mAccuracy;

        /**
         * @brief Maximum number of iterations for the iterative pressure solver.
         */
        int mIterations;
    };
}

int main()
{
    physsim::PhyssimWindow window(
        800,       // width
        600,       // height
        "5_fluid", // title
        std::make_shared<physsim::FluidSimulation>(),
        false // fullscreen
    );

    return window.run();
}
