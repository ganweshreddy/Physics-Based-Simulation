#include <scenes.hpp>

#include <vislab/graphics/bmp_writer.hpp>
#include <vislab/graphics/dielectric_bsdf.hpp>
#include <vislab/graphics/image.hpp>
#include <vislab/graphics/path_radiance_integrator.hpp>
#include <vislab/graphics/sphere.hpp>

int main()
{
    using namespace vislab;

    // image buffer to render into.
    auto image = std::make_shared<Image3d>();
    image->setResolution(128, 128);
    image->setZero();

    // create a scene
    auto scene = pbr::Scenes::create_box_scene(image->getResolution());

    // create sphere 1
    auto sphere1  = std::make_shared<Sphere>();
    sphere1->bsdf = std::make_shared<DielectricBSDF>();
    sphere1->transform.setMatrix(Eigen::Matrix4d::translate(Eigen::Vector3d(0, -0.5, 0)) * Eigen::Matrix4d::scale(Eigen::Vector3d(0.2, 0.2, 0.2)));
    scene->shapes.push_back(sphere1);

    // create sphere 2
    auto sphere2  = std::make_shared<Sphere>();
    sphere2->bsdf = std::make_shared<DielectricBSDF>();
    sphere2->transform.setMatrix(Eigen::Matrix4d::translate(Eigen::Vector3d(-0.3, 0.4, -0.5)) * Eigen::Matrix4d::scale(Eigen::Vector3d(0.3, 0.3, 0.3)));
    scene->shapes.push_back(sphere2);

    // create a light transport integrator and render the image.
    auto integrator = std::make_unique<PathRadianceIntegrator>();
    integrator->setSamplesPerPixel(200);
    integrator->maxDepth      = 4;
    integrator->rouletteDepth = 2;
    integrator->render(scene.get(), image.get());

    // write image to file
    auto writer = std::make_unique<BmpWriter>();
    writer->inputImage.setData(image);
    writer->paramPath.setValue("dielectric.bmp");
    writer->update();

    return 0;
}
