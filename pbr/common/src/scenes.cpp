#include "scenes.hpp"

#include <vislab/graphics/area_light.hpp>
#include <vislab/graphics/const_texture.hpp>
#include <vislab/graphics/diffuse_bsdf.hpp>
#include <vislab/graphics/image.hpp>
#include <vislab/graphics/perspective_camera.hpp>
#include <vislab/graphics/rectangle.hpp>

namespace pbr
{
    std::shared_ptr<vislab::Scene> Scenes::create_box_scene(const Eigen::Vector2i& resolution, bool addLight)
    {
        using namespace vislab;

        // create a scene
        auto scene = std::make_shared<Scene>();

        // create camera
        {
            auto camera = std::make_shared<PerspectiveCamera>();
            camera->setLookAt(Eigen::Vector3d(1, 0, 0));
            camera->setPosition(Eigen::Vector3d(-3.5, 0, 0));
            camera->setUp(Eigen::Vector3d(0, 0, 1));
            camera->setNear(0.01);
            camera->setFar(10);
            camera->setWidth(resolution.x());
            camera->setHeight(resolution.y());
            scene->camera = camera;
        }

        if (addLight)
        {
            // create shape for the area light
            auto rect                   = std::make_shared<Rectangle>();
            rect->bsdf                  = std::make_shared<DiffuseBSDF>(std::make_shared<ConstTexture>(Spectrum(0, 0, 0)));
            Eigen::Matrix4d transform   = Eigen::Matrix4d::Identity();
            transform.block(0, 0, 3, 3) = Eigen::Quaterniond(Eigen::AngleAxisd(EIGEN_PI, Eigen::Vector3d(1, 0, 0))).toRotationMatrix() * Eigen::Matrix3d::Identity() * 0.2;
            transform.block(0, 3, 3, 1) = Eigen::Vector3d(0, 0, 0.9999);
            rect->transform.setMatrix(transform);
            scene->shapes.push_back(rect);

            // create an area light
            auto light      = std::make_shared<AreaLight>();
            light->radiance = Eigen::Vector3d(1, 1, 1) * 10;
            light->shape    = rect;
            rect->light     = light;
            light->transform.setMatrix(Eigen::Matrix4d::translate(Eigen::Vector3d(0.5, 0.6, 0.7)));
            scene->lights.push_back(light);
        }

        {
            // create floor wall
            auto rect  = std::make_shared<Rectangle>();
            rect->bsdf = std::make_shared<DiffuseBSDF>(std::make_shared<ConstTexture>(Spectrum(1, 1, 1)));
            rect->transform.setMatrix(Eigen::Matrix4d::translate(Eigen::Vector3d(0, 0, -1)));
            scene->shapes.push_back(rect);
        }

        {
            // create ceiling wall
            auto rect                   = std::make_shared<Rectangle>();
            rect->bsdf                  = std::make_shared<DiffuseBSDF>(std::make_shared<ConstTexture>(Spectrum(1, 1, 1)));
            Eigen::Matrix4d transform   = Eigen::Matrix4d::Identity();
            transform.block(0, 0, 3, 3) = Eigen::Quaterniond(Eigen::AngleAxisd(EIGEN_PI, Eigen::Vector3d(1, 0, 0))).toRotationMatrix();
            transform.block(0, 3, 3, 1) = Eigen::Vector3d(0, 0, 1);
            rect->transform.setMatrix(transform);
            scene->shapes.push_back(rect);
        }

        {
            // create back wall
            auto rect                   = std::make_shared<Rectangle>();
            rect->bsdf                  = std::make_shared<DiffuseBSDF>(std::make_shared<ConstTexture>(Spectrum(1, 1, 1)));
            Eigen::Matrix4d transform   = Eigen::Matrix4d::Identity();
            transform.block(0, 0, 3, 3) = Eigen::Quaterniond(Eigen::AngleAxisd(-EIGEN_PI / 2, Eigen::Vector3d(0, 1, 0))).toRotationMatrix();
            transform.block(0, 3, 3, 1) = Eigen::Vector3d(1, 0, 0);
            rect->transform.setMatrix(transform);
            scene->shapes.push_back(rect);
        }

        {
            // create right wall
            auto rect                   = std::make_shared<Rectangle>();
            rect->bsdf                  = std::make_shared<DiffuseBSDF>(std::make_shared<ConstTexture>(Spectrum(1, 0.3, 0.3)));
            Eigen::Matrix4d transform   = Eigen::Matrix4d::Identity();
            transform.block(0, 0, 3, 3) = Eigen::Quaterniond(Eigen::AngleAxisd(EIGEN_PI / 2, Eigen::Vector3d(1, 0, 0))).toRotationMatrix();
            transform.block(0, 3, 3, 1) = Eigen::Vector3d(0, 1, 0);
            rect->transform.setMatrix(transform);
            scene->shapes.push_back(rect);
        }

        {
            // create left wall
            auto rect                   = std::make_shared<Rectangle>();
            rect->bsdf                  = std::make_shared<DiffuseBSDF>(std::make_shared<ConstTexture>(Spectrum(0.3, 1.0, 0.3)));
            Eigen::Matrix4d transform   = Eigen::Matrix4d::Identity();
            transform.block(0, 0, 3, 3) = Eigen::Quaterniond(Eigen::AngleAxisd(-EIGEN_PI / 2, Eigen::Vector3d(1, 0, 0))).toRotationMatrix();
            transform.block(0, 3, 3, 1) = Eigen::Vector3d(0, -1, 0);
            rect->transform.setMatrix(transform);
            scene->shapes.push_back(rect);
        }

        return scene;
    }
}
