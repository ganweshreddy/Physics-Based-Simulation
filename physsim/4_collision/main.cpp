#include "collision_detection.hpp"
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
     * @brief Rigid body simulation with collision detection.
     */
    class CollisionSimulation : public Simulation
    {
    public:
        /**
         * @brief Enumeration of integration methods.
         */
        enum EIntegrationMethod
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
            mIntegrationMethod = EIntegrationMethod::ExplicitEuler;
            mBroadPhaseMethod  = EBroadPhaseMethod::SweepAndPrune;
            mNarrowPhaseMethod = ENarrowPhaseMethod::Exhaustive;
            mStepSize          = 3E-3;
            mEpsilon           = 0.5;
            mGravity << 0, -9.81, 0;

            // create a box-shaped mesh
            auto boxMesh       = std::make_shared<Mesh>();
            boxMesh->positions = std::make_shared<Array3f>();
            boxMesh->positions->setValues({ { -1, -1, -1 },
                                            { 1, -1, -1 },
                                            { -1, 1, -1 },
                                            { 1, 1, -1 },
                                            { -1, -1, 1 },
                                            { 1, -1, 1 },
                                            { -1, 1, 1 },
                                            { 1, 1, 1 } });
            // create the index buffer
            boxMesh->indices = std::make_shared<Array3ui>();
            boxMesh->indices->setValues({
                { 0, 1, 4 },
                { 1, 5, 4 },
                { 1, 3, 5 },
                { 3, 7, 5 },
                { 3, 2, 7 },
                { 2, 6, 7 },
                { 2, 0, 6 },
                { 0, 4, 6 },
                { 2, 3, 0 },
                { 3, 1, 0 },
                { 4, 5, 6 },
                { 5, 7, 6 },
            });

            mRigidBodies.resize(NUM_CUBES + 1);
            Eigen::Vector3d c0(204.0 / 255.0, 0, 0);
            Eigen::Vector3d c1(0, 128.0 / 255.0, 204.0 / 255.0);
            for (int i = 0; i < NUM_CUBES; ++i)
            {
                double a    = double(i + 1) / NUM_CUBES;
                auto shape  = std::make_shared<Triangle>();
                shape->mesh = boxMesh;
                shape->bsdf = std::make_shared<DiffuseBSDF>(std::make_shared<ConstTexture>(c0 * a + c1 * (1 - a)));
                scene->shapes.push_back(shape);

                mRigidBodies[i] = std::make_shared<RigidBody>();
                mRigidBodies[i]->setMass(1.);
                mRigidBodies[i]->setInertiaBody(mRigidBodies[i]->mass() * 2.0 / 6.0 * Eigen::Matrix3d::Identity());
                mRigidBodies[i]->setShape(shape);
            }

            // create a box-shaped mesh
            auto groundMesh       = std::make_shared<Mesh>();
            groundMesh->positions = std::make_shared<Array3f>();
            groundMesh->positions->setValues({
                { -1, 0, -1 },
                { 1, 0, -1 },
                { -1, 0, 1 },
                { 1, 0, 1 },
            });
            // create the index buffer
            groundMesh->indices = std::make_shared<Array3ui>();
            groundMesh->indices->setValues({
                { 2, 3, 0 },
                { 3, 1, 0 },
            });

            // create ground plane
            auto ground  = std::make_shared<Triangle>();
            ground->bsdf = std::make_shared<DiffuseBSDF>(std::make_shared<ConstTexture>(Spectrum(1, 1, 1)));
            ground->mesh = groundMesh;
            scene->shapes.push_back(ground);
            mRigidBodies[NUM_CUBES] = std::make_shared<RigidBody>();
            mRigidBodies[NUM_CUBES]->setType(RigidBody::EType::Static);
            mRigidBodies[NUM_CUBES]->setScale(20);
            mRigidBodies[NUM_CUBES]->setPosition({ 0, -5, 0 });
            mRigidBodies[NUM_CUBES]->setShape(ground);
            mRigidBodies[NUM_CUBES]->setMass(100000);
            mRigidBodies[NUM_CUBES]->setInertiaBody(mRigidBodies[NUM_CUBES]->mass() * 2.0 / 6.0 * Eigen::Matrix3d::Identity());

            // create a point light
            auto light       = std::make_shared<PointLight>();
            light->intensity = Eigen::Vector3d(100.0, 100.0, 100.0);
            light->transform.setMatrix(Eigen::Matrix4d::translate(Eigen::Vector3d(-4, 5, 4)));

            // create a camera
            auto camera = std::make_shared<PerspectiveCamera>();
            camera->setLookAt(Eigen::Vector3d(0, 0, 0));
            camera->setPosition(Eigen::Vector3d(-20, 0, -10));
            camera->setUp(Eigen::Vector3d(0, 1, 0));
            camera->setNear(0.01);
            camera->setFar(100);

            // add elements to scene
            scene->camera = camera;
            scene->lights.push_back(light);

            // assign rigid bodies to collision detection solver.
            mCollisionDetection = std::make_unique<CollisionDetection>(mRigidBodies);
        }

        /**
         * @brief Restarts the simulation.
         */
        void restart() override
        {
            double x1 = 5;
            double y1 = 4;
            mRigidBodies[0]->setPosition(Eigen::Vector3d(0, y1, 0));
            mRigidBodies[0]->setRotation(Eigen::Quaterniond(0, -0.3444844, -0.3444844, -0.8733046));
            mRigidBodies[0]->setLinearVelocity(Eigen::Vector3d(std::sin(1.047), std::cos(1.047), 0) * 10);
            mRigidBodies[0]->setAngularVelocity(Eigen::Vector3d::Zero());

            for (size_t i = 1; i < NUM_CUBES; ++i)
            {
                mRigidBodies[i]->setPosition(Eigen::Vector3d(x1, y1 * (i + 1), 0));
                mRigidBodies[i]->setRotation(Eigen::Quaterniond::Identity());
                mRigidBodies[i]->setLinearVelocity(Eigen::Vector3d::Zero());
                mRigidBodies[i]->setAngularVelocity(Eigen::Vector3d::Zero());
            }
        }

        /**
         * @brief Advances the simulation one time step forward.
         * @param elapsedTime Elapsed time in milliseconds during the last frame.
         * @param totalTime Total time in milliseconds since the beginning of the first frame.
         * @param timeStep Time step of the simulation. Restarts when resetting the simulation.
         */
        void advance(double elapsedTime, double totalTime, int64_t timeStep) override
        {
            // set force and torque at beginning of frame to zero
            for (auto body : mRigidBodies)
            {
                body->resetForce();
                body->resetTorque();
            }

            // perform collision detection...
            mCollisionDetection->computeCollisionDetection(
                mBroadPhaseMethod,
                mNarrowPhaseMethod,
                mEpsilon,
                mStepSize);

            // apply gravity
            for (int i = 0; i < NUM_CUBES; ++i)
            {
                auto body = mRigidBodies[i];
                if (body->type() == RigidBody::EType::Static)
                    continue;
                body->applyForceToCenterOfMass(mGravity * body->mass());
            }

            // numerically integrate the bodies
            for (auto body : mRigidBodies)
            {
                switch (mIntegrationMethod)
                {
                case EIntegrationMethod::ExplicitEuler:
                {
                    explicitEuler(*body.get(), mStepSize);
                    break;
                }
                case EIntegrationMethod::SymplecticEuler:
                {
                    symplecticEuler(*body.get(), mStepSize);
                    break;
                }
                case EIntegrationMethod::Implicit:
                {
                    implicitEuler(*body.get(), mStepSize);
                    break;
                }
                }
            }
        }

        /**
         * @brief Adds graphical user interface elements with imgui.
         */
        void gui() override
        {
            ImGui::PushItemWidth(100);

            ImGui::Combo("method", (int*)&mIntegrationMethod, "explicit euler\0symplectic euler\0implicit euler\0\0");
            ImGui::Combo("broad phase", (int*)&mBroadPhaseMethod, "none\0aabb\0sap\0\0");
            ImGui::Combo("narrow phase", (int*)&mNarrowPhaseMethod, "exhaustive\0gjk\0\0");
            
            double stepSizeMin = 1E-3, stepSizeMax = 1E-1;
            ImGui::SliderScalar("dt", ImGuiDataType_Double, &mStepSize, &stepSizeMin, &stepSizeMax);

            double epsilonMin = 1E-3, epsilonMax = 1;
            ImGui::SliderScalar("eps", ImGuiDataType_Double, &mEpsilon, &epsilonMin, &epsilonMax);

            ImGui::PopItemWidth();
        }

    private:
        /**
         * @brief Number of cubes in the scene.
         */
        static const int NUM_CUBES = 5;

        /**
         * @brief Numerical integration method.
         */
        EIntegrationMethod mIntegrationMethod;

        /**
         * @brief Broad phase collision detection method.
         */
        EBroadPhaseMethod mBroadPhaseMethod;

        /**
         * @brief Narrow phase collision detection method.
         */
        ENarrowPhaseMethod mNarrowPhaseMethod;

        /**
         * @brief Rigid bodies in the scene.
         */
        std::vector<std::shared_ptr<RigidBody>> mRigidBodies;

        /**
         * @brief Collision detection solver.
         */
        std::unique_ptr<CollisionDetection> mCollisionDetection;

        /**
         * @brief Integration step size.
         */
        double mStepSize;

        /**
         * @brief Impulse response epsilon.
         */
        double mEpsilon;

        /**
         * @brief Gravitational acceleration.
         */
        Eigen::Vector3d mGravity;
    };
}

int main()
{
    physsim::PhyssimWindow window(
        800,           // width
        600,           // height
        "4_collision", // title
        std::make_shared<physsim::CollisionSimulation>(),
        false // fullscreen
    );

    return window.run();
}
