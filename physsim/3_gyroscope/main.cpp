#include "physsim_window.hpp"
#include "rigid_body.hpp"
#include "rigid_body_integrator.hpp"
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
#include <vislab/graphics/sphere.hpp>
#include <vislab/graphics/triangle.hpp>

using namespace vislab;

namespace physsim
{
    /**
     * @brief Rigid body update with gyroscopic force.
     */
    class GyroscopeSimulation : public Simulation
    {
    public:
        /**
         * @brief Enumeration of integration methods.
         */
        enum EMethod
        {
            ExplicitEuler,
            SymplecticEuler,
            Implicit,
        };

        /**
         * @brief Initializes the scene.
         */
        void init() override
        {
            // initial simulation parameters
            mMethod   = EMethod::ExplicitEuler;
            mStepSize = 1E-2;
            mRigidBody.setMass(1.);
            mRigidBody.setInertiaBody(Eigen::Matrix3d::diag(Eigen::Vector3d(24.1449, 28.436, 118.812)).block(0, 0, 3, 3));

            // create ground plane
            auto rect  = std::make_shared<Rectangle>();
            rect->bsdf = std::make_shared<DiffuseBSDF>(std::make_shared<ConstTexture>(Spectrum(1, 1, 1)));
            Eigen::Matrix4d transform;
            transform << 15, 0, 0, 0,
                0, 0, 15, -3.5,
                0, 15, 0, 0,
                0, 0, 0, 1;
            rect->transform.setMatrix(transform);
            scene->shapes.push_back(rect);

            // create a T-shaped mesh
            auto mesh       = std::make_shared<Mesh>();
            mesh->positions = std::make_shared<Array3f>();
            mesh->positions->setValues({ { -1, -1, -1 },
                                         { 1, -1, -1 },
                                         { -1, 1, -1 },
                                         { 1, 1, -1 },
                                         { -1, -1, 1 },
                                         { 1, -1, 1 },
                                         { -1, 1, 1 },
                                         { 1, 1, 1 },
                                         { -3, -1, -1 },
                                         { 3, -1, -1 },
                                         { -3, 1, -1 },
                                         { 3, 1, -1 },
                                         { -3, -1, 1 },
                                         { 3, -1, 1 },
                                         { -3, 1, 1 },
                                         { 3, 1, 1 },
                                         { -1, -3, -1 },
                                         { 1, -3, -1 },
                                         { -1, -3, 1 },
                                         { 1, -3, 1 } });

            // bring center of mass to (0,0,0)
            auto& mat = mesh->positions->getData();
            auto com  = mesh->positions->getData().rowwise().mean();
            mat.colwise() -= com;

            // create the index buffer
            mesh->indices = std::make_shared<Array3ui>();
            mesh->indices->setValues({
                { 0, 2, 1 },
                { 2, 3, 1 },
                { 9, 11, 13 },
                { 11, 15, 13 },
                { 2, 6, 3 },
                { 6, 7, 3 },
                { 5, 7, 4 },
                { 7, 6, 4 },
                { 12, 14, 8 },
                { 14, 10, 8 },
                { 18, 16, 17 },
                { 18, 17, 19 },
                { 8, 10, 0 },
                { 10, 2, 0 },
                { 3, 9, 1 },
                { 3, 11, 9 },
                { 12, 6, 14 },
                { 12, 4, 6 },
                { 7, 5, 15 },
                { 15, 5, 13 },
                { 14, 2, 10 },
                { 2, 14, 6 },
                { 8, 0, 12 },
                { 4, 12, 0 },
                { 3, 7, 11 },
                { 7, 15, 11 },
                { 5, 1, 9 },
                { 5, 9, 13 },
                { 0, 1, 17 },
                { 0, 17, 16 },
                { 19, 5, 18 },
                { 5, 4, 18 },
                { 0, 16, 4 },
                { 4, 16, 18 },
                { 1, 5, 19 },
                { 17, 1, 19 },
            });

            // create a triangle mesh
            auto shape  = std::make_shared<Triangle>();
            shape->mesh = mesh;
            shape->bsdf = std::make_shared<DiffuseBSDF>(std::make_shared<ConstTexture>(Spectrum(1, 0, 0)));
            mRigidBody.setShape(shape);

            // create a point light
            auto light       = std::make_shared<PointLight>();
            light->intensity = Eigen::Vector3d(100.0, 100.0, 100.0);
            light->transform.setMatrix(Eigen::Matrix4d::translate(Eigen::Vector3d(-4, 5, 4)));

            // create a camera
            auto camera = std::make_shared<PerspectiveCamera>();
            camera->setLookAt(Eigen::Vector3d(0, 0, 0));
            camera->setPosition(Eigen::Vector3d(-10, 0, 0));
            camera->setUp(Eigen::Vector3d(0, 1, 0));
            camera->setNear(0.01);
            camera->setFar(100);

            // add elements to scene
            scene->camera = camera;
            scene->shapes.push_back(shape);
            scene->lights.push_back(light);
        }

        /**
         * @brief Restarts the simulation.
         */
        void restart() override
        {
            mRigidBody.setRotation(Eigen::Quaterniond::Identity());
            mRigidBody.setAngularVelocity(Eigen::Vector3d(0.1, 10, 0));
        }

        /**
         * @brief Advances the simulation one time step forward.
         * @param elapsedTime Elapsed time in milliseconds during the last frame.
         * @param totalTime Total time in milliseconds since the beginning of the first frame.
         * @param timeStep Time step of the simulation. Restarts when resetting the simulation.
         */
        void advance(double elapsedTime, double totalTime, int64_t timeStep) override
        {
            switch (mMethod)
            {
            case EMethod::ExplicitEuler:
            {
                explicitEuler(mRigidBody, mStepSize);
                break;
            }
            case EMethod::SymplecticEuler:
            {
                symplecticEuler(mRigidBody, mStepSize);
                break;
            }
            case EMethod::Implicit:
            {
                implicitEuler(mRigidBody, mStepSize);
                break;
            }
            }
        }

        /**
         * @brief Adds graphical user interface elements with imgui.
         */
        void gui() override
        {
            ImGui::PushItemWidth(100);

            ImGui::Combo("method", (int*)&mMethod, "explicit euler\0symplectic euler\0implicit euler\0\0");
            
            double stepSizeMin = 1E-3, stepSizeMax = 1E-1;
            ImGui::SliderScalar("dt", ImGuiDataType_Double, &mStepSize, &stepSizeMin, &stepSizeMax);

            ImGui::PopItemWidth();
        }

    private:
        /**
         * @brief Numerical integration method.
         */
        EMethod mMethod;

        /**
         * @brief Rigid body state.
         */
        RigidBody mRigidBody;

        /**
         * @brief Integration step size.
         */
        double mStepSize;
    };
}

int main()
{
    physsim::PhyssimWindow window(
        800,           // width
        600,           // height
        "3_gyroscope", // title
        std::make_shared<physsim::GyroscopeSimulation>(),
        false // fullscreen
    );

    return window.run();
}
