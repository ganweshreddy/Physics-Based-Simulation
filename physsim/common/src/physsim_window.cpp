#include "physsim_window.hpp"

#include "imgui_helper.hpp"
#include "renderer.hpp"
#include "simulation.hpp"

#include <vislab/graphics/bmp_writer.hpp>
#include <vislab/graphics/camera.hpp>
#include <vislab/graphics/image.hpp>
#include <vislab/graphics/path_radiance_integrator.hpp>
#include <vislab/graphics/scene.hpp>
#include <vislab/graphics/trackball_interactor.hpp>

#include <backends/imgui_impl_glfw.h>
#include <backends/imgui_impl_opengl3.h>
#include <imgui.h>

namespace physsim
{
    PhyssimWindow::PhyssimWindow(int width, int height, const char* title, std::shared_ptr<Simulation> simulation, bool fullScreen)
        : RenderWindowOpenGL(width, height, title, fullScreen)
        , mSimulation(simulation)
        , mScene(std::make_shared<vislab::Scene>())
        , mRenderer(std::make_shared<Renderer>())
        , mActive(true)
        , mTimeStep(0)
        , mSamplesPerPixel(4)
        , mRecursionDepth(1)
        , mNumFrames(10)
        , mNthTimeStep(10)
    {
    }

    bool PhyssimWindow::init()
    {
        // initialize the renderer and the simulation
        mRenderer->scene   = mScene;
        mSimulation->scene = mScene;
        mSimulation->init();
        mSimulation->restart();

        // connect the interactor with the camera of the scene
        interactor = std::make_shared<vislab::TrackballInteractor>();
        interactor->setCamera(mScene->camera);

        // Setup Dear ImGui context
        IMGUI_CHECKVERSION();
        ImGui::CreateContext();
        ImGuiIO& io = ImGui::GetIO();
        // Setup Platform/Renderer bindings
        ImGui_ImplGlfw_InitForOpenGL(getGlwfWindow(), true);
        ImGui_ImplOpenGL3_Init("#version 330");
        // Setup Dear ImGui style
        ImGui::StyleColorsDark();

        return true;
    }

    void PhyssimWindow::restart()
    {
        mTimeStep = 0;
        mSimulation->restart();
    }

    void PhyssimWindow::advance(double elapsedTime, double totalTime)
    {
        mSimulation->advance(elapsedTime, totalTime, mTimeStep);
        mTimeStep++;
    }

    bool PhyssimWindow::createDevice(vislab::OpenGL* opengl)
    {
        return mRenderer->createDevice(opengl);
    }

    bool PhyssimWindow::createSwapChain(vislab::OpenGL* opengl)
    {
        return mRenderer->createSwapChain(opengl);
    }

    void PhyssimWindow::releaseDevice()
    {
        mRenderer->releaseDevice();
    }

    void PhyssimWindow::releaseSwapChain()
    {
        mRenderer->releaseSwapChain();
    }

    void PhyssimWindow::draw(vislab::OpenGL* opengl, double elapsedTime, double totalTime)
    {
        // advance the simulation one step toward
        if (mActive)
        {
            advance(elapsedTime, totalTime);
        }

        // render the frame with openGL
        mRenderer->draw(opengl);

        // draw UI
        gui();
    }

    void PhyssimWindow::gui()
    {
        // if ImGui wants to process mouse events, disable the interactor
        interactor->active = !ImGui::GetIO().WantCaptureMouse;

        // feed inputs to dear imgui, start new frame
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        // draw UI elements
        ImGui::Begin("PhysSim", 0, ImGuiWindowFlags_AlwaysAutoResize);
        ImguiHelper::toggleButton("active", &mActive);
        ImGui::SameLine();
        if (ImGui::Button("Restart"))
            restart();
        ImGui::Text("time step: %li", mTimeStep);
        ImGui::Separator();

        // simulation UI
        mSimulation->gui();
        ImGui::Separator();

        // path tracing UI
        if (ImGui::CollapsingHeader("Render"))
        {
            ImGui::PushItemWidth(50);
            ImGui::InputInt("spp", &mSamplesPerPixel, 0);
            ImGui::SliderInt("depth", &mRecursionDepth, 1, 8);
            if (ImGui::Button("Path Tracing"))
                pathTrace();
            ImGui::PopItemWidth();
        }
        ImGui::Separator();

        // animation export UI
        if (ImGui::CollapsingHeader("Video"))
        {
            ImGui::PushItemWidth(50);
            ImGui::InputInt("frames", &mNumFrames, 0);
            ImGui::InputInt("nth step", &mNthTimeStep, 0);
            if (ImGui::Button("Export Anim"))
                exportVideo();
            ImGui::PopItemWidth();
        }
        ImGui::Separator();

        // close button
        if (ImGui::Button("Close Window"))
            glfwSetWindowShouldClose(getGlwfWindow(), GL_TRUE);

        ImGui::End();

        // render dear imgui onto screen
        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
    }

    void PhyssimWindow::pathTrace()
    {
        // image buffer to render into.
        auto image = std::make_shared<vislab::Image3d>();
        Eigen::Vector2i resolution(getOpenGL()->getViewport().width, getOpenGL()->getViewport().height);
        image->setResolution(resolution);
        for (int i = 0; i < resolution.prod(); ++i)
            image->setValue(i, Eigen::Vector3d::Ones());

        // create a light transport integrator and render the image.
        auto integrator = std::make_unique<vislab::PathRadianceIntegrator>();
        integrator->setSamplesPerPixel(std::max(0, mSamplesPerPixel));
        integrator->maxDepth      = mRecursionDepth;
        integrator->rouletteDepth = mRecursionDepth;
        integrator->render(mScene.get(), image.get());

        // write image to file
        auto writer = std::make_unique<vislab::BmpWriter>();
        writer->inputImage.setData(image);
        writer->paramPath.setValue("timestep" + std::to_string(mTimeStep) + ".bmp");
        writer->update();
    }

    void PhyssimWindow::exportVideo()
    {
        double dt = 16.666666; // time in milliseconds for 60 fps
        restart();
        for (int i = 0; i < mNumFrames; ++i)
        {
            pathTrace();
            for (int n = 0; n < mNthTimeStep; ++n)
                advance(dt, i * mNthTimeStep * dt);
        }
    }
}
