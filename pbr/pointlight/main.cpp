#include <scenes.hpp>

#include <vislab/graphics/bmp_writer.hpp>
#include <vislab/graphics/image.hpp>
#include <vislab/graphics/path_radiance_integrator.hpp>
#include <vislab/graphics/point_light.hpp>

#include <random>

int main()
{
    using namespace vislab;

    // image buffer to render into.
    auto image = std::make_shared<Image3d>();
    image->setResolution(128, 128);
    image->setZero();

    // create a scene
    auto scene = pbr::Scenes::create_box_scene(image->getResolution(), false);

    // randomly place point lights
    const int numLights = 10;
    std::default_random_engine rng;
    std::uniform_real_distribution<double> rnd(0, 1);
    for (int i = 0; i < numLights; ++i)
    {
        // create a point light
        auto light       = std::make_shared<PointLight>();
        light->intensity = Eigen::Vector3d(rnd(rng), rnd(rng), rnd(rng)) / numLights;
        light->transform.setMatrix(Eigen::Matrix4d::translate(Eigen::Vector3d(rnd(rng) * 1.8 - 0.9, rnd(rng) * 1.8 - 0.9, rnd(rng) * 1.8 - 0.9)));
        scene->lights.push_back(light);
    }

    // create a light transport integrator and render the image.
    auto integrator = std::make_unique<PathRadianceIntegrator>();
    integrator->setSamplesPerPixel(512);
    integrator->maxDepth      = 4;
    integrator->rouletteDepth = 2;
    integrator->render(scene.get(), image.get());

    // write image to file
    auto writer = std::make_unique<BmpWriter>();
    writer->inputImage.setData(image);
    writer->paramPath.setValue("pointlight.bmp");
    writer->update();

    return 0;
}
