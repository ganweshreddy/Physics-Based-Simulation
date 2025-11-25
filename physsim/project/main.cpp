#include "collision_detection.hpp"
#include "physsim_window.hpp"
#include "rigid_body.hpp"
#include "rigid_body_integrator.hpp"
#include "simulation.hpp"

#include <imgui.h>
#include <numbers>
#include <iostream>
#include <vislab/core/array.hpp>
#include <vislab/graphics/const_texture.hpp>
#include <vislab/graphics/diffuse_bsdf.hpp>
#include <vislab/graphics/mesh.hpp>
#include <vislab/graphics/orthographic_camera.hpp>
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
            mNarrowPhaseMethod = ENarrowPhaseMethod::GilbertJohnsonKeerthi;
            mStepSize          = 1E-2;
            mEpsilon           = 0.35;
            mGravity << 0, -9.81, 0;

            // create external sides
            createstaticrect(10, 200, 2, { -55 + 1.75 - .25, 100, 0 }, true); // right
            createstaticrect(10, 200, 2, { 55 - 1.75 + .25, 100, 0 }, true);  // left
            createstaticrect(117.5, 10, 2, { 0, -5, 0 }, true);               // bottom
            createstaticrect(116.5, 210, 2, { 0, 210 / 2 - 10, 2 }, false);   // front
            createstaticrect(116.5, 210, 2, { 0, 210 / 2 - 10, -2 }, false);  // back
            createstaticrect(5, 60, 2, { -25, 220, 0 }, true, 10);            // top right
            createstaticrect(5, 60, 2, { 25, 220, 0 }, true, -10);            // top left

            // create deviders
            for (size_t i = 0; i < 11; ++i)
                createstaticrect(4.5, 50, 2, { 42.5 - 8.5 * i, 25, 0 }, true);

            // create the pegs
            for (size_t i = 0; i < 7; ++i)
            {
                createstaticpattern({ -44.0, 58.0 + 20 * i, -1.0 }, 9, 11);
                createstaticpattern({ -39.75, 66.5 + 20 * i, -1.0 }, 8, 11);
            }

            // create the dynamic balls

            double radius = 2;

            ballstart = createcylinder(radius, { -2, 190, -1 }, physsim::RigidBody::EType::Dynamic, 1);
            /*
            int xdimballs = sqrt(numballs);
            int ydimballs = numballs / xdimballs;
            for (size_t i = 0; i < xdimballs; ++i)
                for (size_t j = 0; j < ydimballs; j++)
                    ballend = createcylinder(radius, { 20 + i * 4.5, 190 + j * 4.5, -1 }, physsim::RigidBody::EType::Dynamic, 1);
            */
            for (size_t i = 1; i < numballs; i++)
            {
                ballend = createcylinder(radius, { -2, 190 + i * 4.5, -1 }, physsim::RigidBody::EType::Dynamic, 1);
            }

            // create a point light
            auto light       = std::make_shared<PointLight>();
            light->intensity = Eigen::Vector3d(1e5, 1e5, 1e5);
            light->transform.setMatrix(Eigen::Matrix4d::translate(Eigen::Vector3d(0, 95, 300)));

            // create a camera
            auto camera = std::make_shared<PerspectiveCamera>();
            camera->setLookAt(Eigen::Vector3d(0, 95, -13));
            camera->setPosition(Eigen::Vector3d(0, 100, 400));
            camera->setUp(Eigen::Vector3d(0, 1, 0));
            camera->setNear(50);
            camera->setFar(500);

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
            Eigen::Vector3d startloc = { -20, 230, -1 };
            for (size_t i = 0; i < numballs; ++i)
            {
                mRigidBodies[i + ballstart]->setPosition(startloc + 4.5 * Eigen::Vector3d(i / pattern, i % pattern, 0));
                mRigidBodies[i + ballstart]->setRotation(Eigen::Quaterniond::Identity());
                mRigidBodies[i + ballstart]->setLinearVelocity(Eigen::Vector3d::Zero());
                mRigidBodies[i + ballstart]->setAngularVelocity(Eigen::Vector3d::Zero());
                mRigidBodies[i + ballstart]->resetForce();
                mRigidBodies[i + ballstart]->resetTorque();
                mRigidBodies[i + ballstart]->setType(RigidBody::EType::Dynamic);
                auto shape  = std::const_pointer_cast<Shape>(mRigidBodies[i + ballstart]->shape());
                shape->bsdf = std::make_shared<DiffuseBSDF>(std::make_shared<ConstTexture>(Eigen::Vector3d(1, 1, 1)));
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
            bool timedo = !true;
            #if timedo
            auto start = std::chrono::system_clock::now();
            #endif
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
            for (auto& body : mRigidBodies)
            {
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

            for (auto& body : mRigidBodies)
            {
                if (body->type() == RigidBody::EType::Static)
                {
                    continue;
                    body->resetForce();
                    body->resetTorque();
                }
                if (body->position()(1) < 50.0)
                {
                    mStepSize = 1e-2;
                    if (body->veldiff() < 5e-2)
                    {
                        body->setType(RigidBody::EType::Static);
                        body->resetForce();
                        body->resetTorque();
                        auto shape  = std::const_pointer_cast<Shape>(body->shape());
                        shape->bsdf = std::make_shared<DiffuseBSDF>(std::make_shared<ConstTexture>(Eigen::Vector3d(1, 0, 0)));
                    }
                }
                // check the limit of val to set static
            }
            #if timedo
            auto end     = std::chrono::system_clock::now();
            auto elapsed = end - start;
            std::cout << elapsed.count() << '\n';
            #endif
        }
        /*
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
        void createstaticpattern(Eigen::Vector3d startloc, size_t num, float gap = 8.5)
        {
            auto val{ 0 };
            for (size_t i = 0; i < num; i++)
                val = createcylinder(4.5 / 2, startloc + Eigen::Vector3d(gap * i, 0, 0), physsim::RigidBody::EType::Static);
        }

        size_t createstaticrect(float x, float y, float z, Eigen::Vector3d loc, bool toscene, double angle = 0)
        {
            auto sides       = createBoxMesh(x / 2, y / 2, z / 2);
            auto sidesshape  = std::make_shared<Triangle>();
            sidesshape->bsdf = std::make_shared<DiffuseBSDF>(std::make_shared<ConstTexture>(Eigen::Vector3d(1, 1, 1)));
            sidesshape->mesh = sides;
            if (toscene)
                scene->shapes.push_back(sidesshape);
            mRigidBodies.push_back(std::make_shared<RigidBody>());
            mRigidBodies.back()->setType(RigidBody::EType::Static);
            mRigidBodies.back()->setScale(1);
            mRigidBodies.back()->setPosition(loc);
            mRigidBodies.back()->setRotation(Eigen::Quaterniond(Eigen::AngleAxisd(angle, Eigen::Vector3d::UnitZ())));
            mRigidBodies.back()->setShape(sidesshape);
            mRigidBodies.back()->setMass(1e10);
            mRigidBodies.back()->setInertiaBody(mRigidBodies[0]->mass() * 2.0 / 6.0 * Eigen::Matrix3d::Identity());
            return mRigidBodies.size() - 1;
        }

        size_t createcylinder(float radius, Eigen::Vector3d loc, const physsim::RigidBody::EType type, double mass = 1e10)
        {
            auto sides  = createCylindermesh(radius);
            auto shape  = std::make_shared<Triangle>();
            auto a      = ballstart / mRigidBodies.size();
            shape->bsdf = std::make_shared<DiffuseBSDF>(std::make_shared<ConstTexture>(Eigen::Vector3d(1, 1, 1)));
            shape->mesh = sides;
            scene->shapes.push_back(shape);
            mRigidBodies.push_back(std::make_shared<RigidBody>());
            mRigidBodies.back()->setType(type);
            mRigidBodies.back()->setScale(1);
            mRigidBodies.back()->setPosition(loc);
            mRigidBodies.back()->setShape(shape);
            mRigidBodies.back()->setMass(mass);
            mRigidBodies.back()->setInertiaBody(mRigidBodies[0]->mass() * 2.0 / 6.0 * Eigen::Matrix3d::Identity());
            return mRigidBodies.size() - 1;
        }

        std::shared_ptr<Mesh> createBoxMesh(float halfdepth = 1, float halfht = 1, float halthk = 1)
        {
            // create a box-shaped mesh
            auto boxMesh       = std::make_shared<Mesh>();
            boxMesh->positions = std::make_shared<Array3f>();

            int neg1{}, neg2{}, neg3{};
            for (size_t i = 0; i < 8; i++)
            {
                neg1 = (i % 2 != 0) ? -1 : +1;
                neg2 = (i % 4 < 2) ? -1 : +1;
                neg3 = (i < 4) ? -1 : +1;
                boxMesh->positions->append({ neg1 * halfdepth, neg2 * halfht, neg3 * halthk });
            }
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
            return boxMesh;
        }

        std::shared_ptr<Mesh> createCylindermesh(float radius)
        {
            auto cylMesh           = std::make_shared<Mesh>();
            cylMesh->positions     = std::make_shared<Array3f>();
            const int num          = 40;
            float h1               = 0;
            float h2               = 2;
            float circley[num + 2] = {};
            float circlez[num + 2] = {};
            // circle 1
            cylMesh->positions->append({ 0, 0, 0 });
            for (size_t i = 0; i <= num; i++)
            {
                auto angle     = i * (3.141592653 * 2) / num;
                circley[i + 1] = radius * cos(angle);
                circlez[i + 1] = radius * sin(angle);

                cylMesh->positions->append({ circley[i + 1], circlez[i + 1], 0 });
            }

            // circle2

            cylMesh->positions->append({ 0, 0, 2 });

            for (size_t i = 0; i <= num; i++)
            {

                cylMesh->positions->append({ circley[i + 1], circlez[i + 1], 2 });
            }

            // nodes of circle

            unsigned int circle1nodes[num + 2];
            unsigned int circle2nodes[num + 2];

            for (unsigned int i = 0; i < num + 2; i++)
            {
                circle1nodes[i] = i;
                circle2nodes[i] = i + num + 2;
            }

            // circle1

            cylMesh->indices = std::make_shared<Array3ui>();
            for (int i = 0; i < num; i++)
            {
                cylMesh->indices->append({ circle1nodes[0], circle1nodes[i + 1], circle1nodes[i + 2] });
            }

            // circle2

            for (int i = 0; i < num; i++)
            {
                cylMesh->indices->append({ circle2nodes[0], circle2nodes[i + 1], circle2nodes[i + 2] });
            }

            // circle1to2

            for (unsigned int i = 1; i <= num; i++)
            {
                cylMesh->indices->append({ circle1nodes[i], circle1nodes[i + 1], circle2nodes[i] });
            }

            // circle2to1

            for (unsigned int i = 1; i <= num; i++)
            {
                cylMesh->indices->append({ circle1nodes[i + 1], circle2nodes[i], circle2nodes[i + 1] });
            }

            return cylMesh;
        }

        size_t ballstart = 0;
        size_t ballend   = 0;
        size_t pattern   = 9;
        size_t numballs  = std::pow(pattern, 2);

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
        800,            // width
        600,            // height
        "Galton Board", // title
        std::make_shared<physsim::CollisionSimulation>(),
        false // fullscreen
    );

    return window.run();
}
