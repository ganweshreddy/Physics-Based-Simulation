#include "renderer.hpp"

#include <vislab/core/array.hpp>
#include <vislab/field/regular_field.hpp>
#include <vislab/graphics/colormap_texture.hpp>
#include <vislab/graphics/const_texture.hpp>
#include <vislab/graphics/diffuse_bsdf.hpp>
#include <vislab/graphics/mesh.hpp>
#include <vislab/graphics/point_light.hpp>
#include <vislab/graphics/projective_camera.hpp>
#include <vislab/graphics/rectangle.hpp>
#include <vislab/graphics/scene.hpp>
#include <vislab/graphics/sphere.hpp>
#include <vislab/graphics/triangle.hpp>
#include <vislab/opengl/opengl.hpp>
#include <vislab/opengl/shader.hpp>

namespace physsim
{
    bool Renderer::createDevice(vislab::OpenGL* opengl)
    {
        if (!mParams.createDevice())
            return false;

        // create a VBO and VAO containing one point
        GLuint sphereVBO;
        glGenBuffers(1, &sphereVBO);
        glGenVertexArrays(1, &mSphereVAO);
        float spherePos[] = { 0.0f, 0.0f, 0.0f };
        glBindVertexArray(mSphereVAO);
        glEnableVertexAttribArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, sphereVBO);
        glBufferData(GL_ARRAY_BUFFER, 3 * sizeof(float), &spherePos, GL_DYNAMIC_DRAW);
        glVertexAttribPointer(0, 3, GL_FLOAT, false, 0, 0);
        glBindVertexArray(0);

        // create a VBO and VAO containing one point
        GLuint rectangleVBOPos, rectangleVBOTexCoord;
        glGenBuffers(1, &rectangleVBOPos);
        glGenBuffers(1, &rectangleVBOTexCoord);
        glGenVertexArrays(1, &mRectangleVAO);
        float rectanglePos[]      = { -1.f, -1.f, 0.f, 1.f, -1.f, 0.f, -1.f, 1.f, 0.f, 1.f, 1.f, 0.f };
        float rectangleTexCoord[] = { 0, 0, 1, 0, 0, 1, 1, 1 };
        glBindVertexArray(mRectangleVAO);
        glBindBuffer(GL_ARRAY_BUFFER, rectangleVBOPos);
        glEnableVertexAttribArray(0);
        glBufferData(GL_ARRAY_BUFFER, 12 * sizeof(float), &rectanglePos, GL_DYNAMIC_DRAW);
        glVertexAttribPointer(0, 3, GL_FLOAT, false, 0, 0);
        glBindBuffer(GL_ARRAY_BUFFER, rectangleVBOTexCoord);
        glEnableVertexAttribArray(1);
        glBufferData(GL_ARRAY_BUFFER, 8 * sizeof(float), &rectangleTexCoord, GL_DYNAMIC_DRAW);
        glVertexAttribPointer(1, 2, GL_FLOAT, false, 0, 0);
        glBindVertexArray(0);

        std::string cb =
            "layout(std140) uniform cbparam\n"
            "{\n"
            "	mat4 worldMatrix;\n"
            "	mat4 normalMatrix;\n"
            "	mat4 viewMatrix;\n"
            "	mat4 projMatrix;\n"
            "	mat4 invMVPMatrix;\n"
            "	vec4 eyePosition;\n"
            "	vec4 diffuse;\n"
            "	vec4 lightPos[8];\n"
            "	vec4 lightIntensity[8];\n"
            "};\n";

        std::string sphere_vs =
            "#version 330 core\n"
            "layout (location = 0) in vec3 Position;\n" +
            cb +
            "void main()\n"
            "{\n"
            "   gl_Position = worldMatrix * vec4(Position, 1.0);\n"
            "}\0";

        std::string sphere_gs =
            "#version 330 core\n"
            "layout(points) in;\n"
            "layout(triangle_strip, max_vertices = 4) out;\n"
            "out vec2 gs_out_TexCoord;\n"
            "out vec4 gs_out_Position;\n"
            "out vec3 gs_out_Center;\n" +
            cb +
            "void main()\n"
            "{\n"
            "   gs_out_Center = gl_in[0].gl_Position.xyz;\n"
            "\n"
            "   vec3 eye = normalize(eyePosition.xyz - gs_out_Center);\n"
            "   vec4 viewPos = viewMatrix * vec4(gs_out_Center, 1);\n"
            "   viewPos /= viewPos.w;\n"
            "\n"
            "   gl_Position = projMatrix * (viewPos + vec4(-1, -1, 0, 0) * sqrt(2));\n"
            "   gs_out_Position = gl_Position;\n"
            "   gs_out_TexCoord = vec2(0,0);\n"
            "   EmitVertex();\n"
            "\n"
            "   gl_Position = projMatrix * (viewPos + vec4(1, -1, 0, 0) * sqrt(2));\n"
            "   gs_out_Position = gl_Position;\n"
            "   gs_out_TexCoord = vec2(1,0);\n"
            "   EmitVertex();\n"
            "\n"
            "   gl_Position = projMatrix * (viewPos + vec4(-1, 1, 0, 0) * sqrt(2));\n"
            "   gs_out_Position = gl_Position;\n"
            "   gs_out_TexCoord = vec2(0,1);\n"
            "   EmitVertex();\n"
            "\n"
            "   gl_Position = projMatrix * (viewPos + vec4(1, 1, 0, 0) * sqrt(2));\n"
            "   gs_out_Position = gl_Position;\n"
            "   gs_out_TexCoord = vec2(1,1);\n"
            "   EmitVertex();\n"
            "   EndPrimitive();\n"
            "}\0";

        std::string intersectRaySphere =
            "bool intersectRaySphere(vec3 o, vec3 l, out vec2 d)\n"
            "{\n"
            "   float dlo = dot(l, o);\n"
            "   float disc = pow(dlo, 2) - dot(o, o) + 1;\n"
            "   if (disc < 0) return false;\n"
            "   float sqrtdisc = sqrt(disc);\n"
            "   d = vec2(-dlo - sqrtdisc, -dlo + sqrtdisc);\n"
            "   return true;\n"
            "}\n";

        std::string sphere_ps =
            "#version 330 core\n"
            "in vec2 gs_out_TexCoord;\n"
            "in vec4 gs_out_Position;\n"
            "in vec3 gs_out_Center;\n"
            "out vec4 FragColor;\n" +
            cb + intersectRaySphere +
            "void main()\n"
            "{\n"
            "   vec4 npc = gs_out_Position;\n"
            "   npc /= npc.w;\n"
            "   vec4 npc0 = vec4(npc.xy, -1, 1);\n"
            "   vec4 npc1 = vec4(npc.xy, 1, 1);\n"
            "   vec4 local0 = invMVPMatrix * npc0;\n"
            "   vec4 local1 = invMVPMatrix * npc1;\n"
            "   vec3 origin = local0.xyz / local0.w;\n"
            "   vec3 dir = normalize(local1.xyz / local1.w - origin);\n"
            "   vec2 d;\n"
            "   if (!intersectRaySphere(origin, dir, d))\n"
            "      discard;\n"
            "   float t = d.x < 0 ? d.y : d.x;\n"
            "   vec4 worldx = worldMatrix * vec4(origin + dir * t, 1);\n"
            "   vec4 npcx = projMatrix * viewMatrix * worldx;\n"
            "   gl_FragDepth = npcx.z / npcx.w * 0.5 + 0.5;\n"
            "   // compute shading\n"
            "   vec3 N = normalize(worldx.xyz - gs_out_Center);\n"
            "   vec3 color = vec3(0,0,0);\n"
            "   for (int i = 0; i < 8; ++i) {\n"
            "      vec3 L = lightPos[i].xyz - worldx.xyz;\n"
            "      float dist2 = dot(L, L);\n"
            "      float cos_theta_o = max(0, dot(N, normalize(L)));\n"
            "      color += diffuse.rgb / 3.1415926535 * lightIntensity[i].rgb * cos_theta_o / dist2;\n"
            "   }\n"
            "   FragColor = vec4(color, diffuse.a);\n"
            "}\0";

        mSphereShader = std::make_shared<vislab::ShaderOpenGL>();
        mSphereShader->setVertexShaderSource(sphere_vs);
        mSphereShader->setGeometryShaderSource(sphere_gs);
        mSphereShader->setPixelShaderSource(sphere_ps);
        if (!mSphereShader->createDevice())
            return false;

        std::string rectangle_vs =
            "#version 330 core\n"
            "layout (location = 0) in vec3 Position;\n"
            "out vec4 vs_out_Position;\n" +
            cb +
            "void main()\n"
            "{\n"
            "   vs_out_Position = worldMatrix * vec4(Position, 1.0);\n"
            "   gl_Position = projMatrix * viewMatrix * vs_out_Position;\n"
            "}\0";

        std::string rectangle_ps =
            "#version 330 core\n"
            "in vec4 vs_out_Position;\n"
            "out vec4 FragColor;\n" +
            cb +
            "void main()\n"
            "{\n"
            "   vec3 N = normalize((worldMatrix * vec4(0,0,1,0)).xyz);\n"
            "   vec3 color = vec3(0,0,0);\n"
            "   for (int i = 0; i < 8; ++i) {\n"
            "      vec3 L = lightPos[i].xyz - vs_out_Position.xyz;\n"
            "      float dist2 = dot(L, L);\n"
            "      float cos_theta_o = max(0, dot(N, normalize(L)));\n"
            "      color += diffuse.rgb / 3.1415926535 * lightIntensity[i].rgb * cos_theta_o / dist2;\n"
            "   }\n"
            "   FragColor = vec4(color, diffuse.a);\n"
            "}\0";

        mRectangleShader = std::make_shared<vislab::ShaderOpenGL>();
        mRectangleShader->setVertexShaderSource(rectangle_vs);
        mRectangleShader->setPixelShaderSource(rectangle_ps);
        if (!mRectangleShader->createDevice())
            return false;

        std::string textured_rectangle_vs =
            "#version 330 core\n"
            "layout (location = 0) in vec3 Position;\n"
            "layout (location = 1) in vec2 TexCoord;\n"
            "out vec4 vs_out_Position;\n"
            "out vec2 vs_out_TexCoord;\n" +
            cb +
            "void main()\n"
            "{\n"
            "   vs_out_Position = worldMatrix * vec4(Position, 1.0);\n"
            "   vs_out_TexCoord = TexCoord;\n"
            "   gl_Position = projMatrix * viewMatrix * vs_out_Position;\n"
            "}\0";

        std::string textured_rectangle_ps =
            "#version 330 core\n"
            "in vec4 vs_out_Position;\n"
            "in vec2 vs_out_TexCoord;\n"
            "uniform sampler2D in_Texture;\n"
            "out vec4 FragColor;\n" +
            cb +
            "void main()\n"
            "{\n"
            "   FragColor = texture(in_Texture, vs_out_TexCoord);\n"
            "}\0";

        mTexturedRectangleShader = std::make_shared<vislab::ShaderOpenGL>();
        mTexturedRectangleShader->setVertexShaderSource(textured_rectangle_vs);
        mTexturedRectangleShader->setPixelShaderSource(textured_rectangle_ps);
        if (!mTexturedRectangleShader->createDevice())
            return false;

        std::string flat_triangle_vs =
            "#version 330 core\n"
            "layout (location = 0) in vec3 Position;\n" +
            cb +
            "void main()\n"
            "{\n"
            "   gl_Position = worldMatrix * vec4(Position, 1.0);\n"
            "}\0";

        std::string flat_triangle_gs =
            "#version 330 core\n"
            "layout(triangles) in;\n"
            "layout(triangle_strip, max_vertices = 3) out;\n"
            "out vec3 gs_out_Position;\n"
            "out vec3 gs_out_Normal;\n" +
            cb +
            "void main()\n"
            "{\n"
            "   gs_out_Normal = normalize(cross(gl_in[1].gl_Position.xyz - gl_in[0].gl_Position.xyz, gl_in[2].gl_Position.xyz - gl_in[0].gl_Position.xyz));\n"
            "   for (int i=0; i<3; ++i) {\n"
            "      gs_out_Position = gl_in[i].gl_Position.xyz;\n"
            "      gl_Position = projMatrix * viewMatrix * vec4(gs_out_Position, 1);\n"
            "      EmitVertex();\n"
            "   }\n"
            "   EndPrimitive();\n"
            "}\0";

        std::string flat_triangle_ps =
            "#version 330 core\n"
            "in vec3 gs_out_Position;\n"
            "in vec3 gs_out_Normal;\n"
            "out vec4 FragColor;\n" +
            cb +
            "void main()\n"
            "{\n"
            "   vec3 N = normalize(gs_out_Normal);\n"
            "   vec3 color = vec3(0,0,0);\n"
            "   for (int i = 0; i < 8; ++i) {\n"
            "      vec3 L = lightPos[i].xyz - gs_out_Position;\n"
            "      float dist2 = dot(L, L);\n"
            "      float cos_theta_o = max(0, dot(N, normalize(L)));\n"
            "      color += diffuse.rgb / 3.1415926535 * lightIntensity[i].rgb * cos_theta_o / dist2;\n"
            "   }\n"
            "   FragColor = vec4(color, diffuse.a);\n"
            "}\0";

        mFlatTriangleShader = std::make_shared<vislab::ShaderOpenGL>();
        mFlatTriangleShader->setVertexShaderSource(flat_triangle_vs);
        mFlatTriangleShader->setGeometryShaderSource(flat_triangle_gs);
        mFlatTriangleShader->setPixelShaderSource(flat_triangle_ps);
        if (!mFlatTriangleShader->createDevice())
            return false;

        std::string triangle_vs =
            "#version 330 core\n"
            "layout (location = 0) in vec3 Position;\n"
            "layout (location = 1) in vec3 Normal;\n"
            "out vec3 gs_out_Position;\n"
            "out vec3 gs_out_Normal;\n" +
            cb +
            "void main()\n"
            "{\n"
            "   gs_out_Position = (worldMatrix * vec4(Position, 1.0)).xyz;\n"
            "   gl_Position = projMatrix * viewMatrix * vec4(gs_out_Position, 1);"
            "   gs_out_Normal = (normalMatrix * vec4(Normal, 0.0)).xyz;\n"
            "}\0";

        mTriangleShader = std::make_shared<vislab::ShaderOpenGL>();
        mTriangleShader->setVertexShaderSource(triangle_vs);
        mTriangleShader->setPixelShaderSource(flat_triangle_ps);
        if (!mTriangleShader->createDevice())
            return false;

        return true;
    }

    bool Renderer::createSwapChain(vislab::OpenGL* opengl)
    {
        return true;
    }

    void Renderer::releaseDevice()
    {
        mParams.releaseDevice();
        mSphereShader->releaseDevice();
        mRectangleShader->releaseDevice();
        mFlatTriangleShader->releaseDevice();
        mTriangleShader->releaseDevice();
        glDeleteVertexArrays(1, &mSphereVAO);
        glDeleteVertexArrays(1, &mRectangleVAO);
        for (auto& t : mTriangleMeshes)
        {
            glDeleteBuffers(1, &t.second.vbo_positions);
            glDeleteBuffers(1, &t.second.vbo_normals);
            glDeleteBuffers(1, &t.second.ibo);
            glDeleteVertexArrays(1, &t.second.vao);
        }
        for (auto& t : mTextures)
        {
            unsigned int tex = t.second;
            glDeleteTextures(1, &tex);
        }
    }

    void Renderer::releaseSwapChain()
    {
    }

    void Renderer::draw(vislab::OpenGL* opengl)
    {
        // real-time renderer only supports projective camera
        auto projectiveCamera = std::dynamic_pointer_cast<vislab::ProjectiveCamera>(scene->camera);
        if (!projectiveCamera)
            return;

        // set camera viewport
        projectiveCamera->setWidth(opengl->getViewport().width);
        projectiveCamera->setHeight(opengl->getViewport().height);

        // update the light sources in the uniform buffer
        int ilight = 0;
        for (auto light : scene->lights)
        {
            auto pointLight = std::dynamic_pointer_cast<vislab::PointLight>(light);
            if (pointLight)
            {
                mParams.data.lightIntensity[ilight] = Eigen::Vector4f(
                    static_cast<float>(pointLight->intensity.x()),
                    static_cast<float>(pointLight->intensity.y()),
                    static_cast<float>(pointLight->intensity.z()),
                    1.);
                mParams.data.lightPos[ilight] = Eigen::Vector4f(
                    static_cast<float>(pointLight->transform.getMatrix().translation().x()),
                    static_cast<float>(pointLight->transform.getMatrix().translation().y()),
                    static_cast<float>(pointLight->transform.getMatrix().translation().z()),
                    1.);
                ilight++;
                if (ilight == 8)
                    break;
            }
        }
        for (int i = ilight; i < 8; ++i)
        {
            mParams.data.lightIntensity[i] = Eigen::Vector4f::Zero();
            mParams.data.lightPos[i]       = Eigen::Vector4f::Zero();
        }

        // clear the backbuffer and depth buffer
        Eigen::Vector4f clearColor(1, 1, 1, 1);
        glClearColor(clearColor.x(), clearColor.y(), clearColor.z(), clearColor.w());
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glEnable(GL_FRAMEBUFFER_SRGB);

        // setup depth writing and test
        glEnable(GL_DEPTH_TEST);
        glDepthMask(GL_TRUE);
        glDepthFunc(GL_LESS);

        // iterate all scene elements
        for (auto shape : scene->shapes)
        {
            // set model view projection matrix for the current object
            mParams.data.worldMatrix  = shape->transform.getMatrix().cast<float>();
            mParams.data.normalMatrix = Eigen::Matrix4f::Identity();
            mParams.data.normalMatrix.block(0, 0, 3, 3) = shape->transform.getMatrixInverse().block(0, 0, 3, 3).transpose().cast<float>();
            mParams.data.viewMatrix   = projectiveCamera->getView().cast<float>();
            mParams.data.projMatrix   = projectiveCamera->getProj().cast<float>();
            mParams.data.eyePosition  = Eigen::Vector4d(projectiveCamera->getPosition().x(), projectiveCamera->getPosition().y(), projectiveCamera->getPosition().z(), 1.).cast<float>();
            mParams.data.invMVPMatrix = (mParams.data.projMatrix * mParams.data.viewMatrix * mParams.data.worldMatrix).inverse();

            // set bsdf parameters
            auto diffuse = std::dynamic_pointer_cast<vislab::DiffuseBSDF>(shape->bsdf);
            if (diffuse)
            {
                auto constTex = std::dynamic_pointer_cast<vislab::ConstTexture>(diffuse->reflectance);
                if (constTex)
                    mParams.data.diffuse = Eigen::Vector4f(
                        static_cast<float>(constTex->color.x()),
                        static_cast<float>(constTex->color.y()),
                        static_cast<float>(constTex->color.z()),
                        1.f);
                else
                    mParams.data.diffuse = Eigen::Vector4f::Ones();
            }
            else
                mParams.data.diffuse = Eigen::Vector4f::Ones();

            // copy parameters to GPU
            mParams.updateBuffer();

            // is a sphere shape?
            auto sphere = std::dynamic_pointer_cast<vislab::Sphere>(shape);
            if (sphere)
            {
                mSphereShader->bind();
                mParams.bind(mSphereShader->getProgram(), "cbparam", 0);
                glBindVertexArray(mSphereVAO);
                glDrawArrays(GL_POINTS, 0, 1);
                glBindVertexArray(0);
                glUseProgram(0);
            }

            // is a rectangle shape?
            auto rectangle = std::dynamic_pointer_cast<vislab::Rectangle>(shape);
            if (rectangle)
            {
                auto colorTex = std::dynamic_pointer_cast<vislab::ColormapTexture>(diffuse->reflectance);
                if (colorTex)
                {
                    unsigned texres = getTextureResources(colorTex.get());
                    glActiveTexture(GL_TEXTURE0);
                    glBindTexture(GL_TEXTURE_2D, texres);
                    mTexturedRectangleShader->bind();
                    mParams.bind(mTexturedRectangleShader->getProgram(), "cbparam", 0);
                    glUniform1i(glGetUniformLocation(mTexturedRectangleShader->getProgram(), "in_Texture"), 0);
                }
                else
                {
                    mRectangleShader->bind();
                    mParams.bind(mRectangleShader->getProgram(), "cbparam", 0);
                }
                glBindVertexArray(mRectangleVAO);
                glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
                glBindVertexArray(0);
                glUseProgram(0);
            }

            // is a triangle shape?
            auto triangle = std::dynamic_pointer_cast<vislab::Triangle>(shape);
            if (triangle)
            {
                if (triangle->mesh->normals)
                {
                    mTriangleShader->bind();
                    mParams.bind(mTriangleShader->getProgram(), "cbparam", 0);
                }
                else
                {
                    mFlatTriangleShader->bind();
                    mParams.bind(mFlatTriangleShader->getProgram(), "cbparam", 0);
                }
                auto bos = getMeshResources(triangle->mesh.get());
                glBindVertexArray(bos.vao);
                glDrawElements(GL_TRIANGLES, triangle->mesh->indices->getSize() * 3, GL_UNSIGNED_INT, 0);
                glBindVertexArray(0);
                glUseProgram(0);
            }
        }
    }

    void Renderer::clearMeshResources(vislab::Mesh* mesh)
    {
        auto it = mTriangleMeshes.find(mesh);
        if (it != mTriangleMeshes.end())
        {
            glDeleteBuffers(1, &it->second.vbo_positions);
            glDeleteBuffers(1, &it->second.vbo_normals);
            glDeleteBuffers(1, &it->second.ibo);
            glDeleteVertexArrays(1, &it->second.vao);
            mTriangleMeshes.erase(mesh);
        }
    }

    void Renderer::clearTextureResources(vislab::ColormapTexture* texture)
    {
        auto it = mTextures.find(texture);
        if (it != mTextures.end())
        {
            unsigned int tex = it->second;
            glDeleteTextures(1, &tex);
            mTextures.erase(texture);
        }
    }

    Renderer::MeshResource Renderer::getMeshResources(vislab::Mesh* mesh)
    {
        auto it = mTriangleMeshes.find(mesh);
        if (it != mTriangleMeshes.end())
            return it->second;

        if (!mesh->positions || !mesh->indices)
            throw std::logic_error("Vertex buffer and index buffer needed!");

        // create VBO for positions
        GLuint triangleVBO_positions;
        glGenBuffers(1, &triangleVBO_positions);
        glBindBuffer(GL_ARRAY_BUFFER, triangleVBO_positions);
        glBufferData(GL_ARRAY_BUFFER, mesh->positions->getSizeInBytes(), mesh->positions->getData().data(), GL_DYNAMIC_DRAW);

        // create VBO for normals
        GLuint triangleVBO_normals = -1;
        if (mesh->normals)
        {
            glGenBuffers(1, &triangleVBO_normals);
            glBindBuffer(GL_ARRAY_BUFFER, triangleVBO_normals);
            glBufferData(GL_ARRAY_BUFFER, mesh->normals->getSizeInBytes(), mesh->normals->getData().data(), GL_DYNAMIC_DRAW);
        }

        // create IBO
        GLuint triangleIBO;
        glGenBuffers(1, &triangleIBO);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, triangleIBO);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, mesh->indices->getSizeInBytes(), mesh->indices->getData().data(), GL_DYNAMIC_DRAW);

        // create VAO
        GLuint triangleVAO;
        glGenVertexArrays(1, &triangleVAO);
        glBindVertexArray(triangleVAO);
        glEnableVertexAttribArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, triangleVBO_positions);
        glVertexAttribPointer(0, 3, GL_FLOAT, false, 0, 0);
        if (mesh->normals)
        {
            glEnableVertexAttribArray(1);
            glBindBuffer(GL_ARRAY_BUFFER, triangleVBO_normals);
            glVertexAttribPointer(1, 3, GL_FLOAT, false, 0, 0);
        }
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, triangleIBO);
        glBindVertexArray(0);

        // memorize the tuple
        MeshResource resource = {
            triangleVBO_positions,
            triangleVBO_normals,
            triangleIBO,
            triangleVAO
        };
        mTriangleMeshes.insert(std::make_pair(mesh, resource));

        // hook an event handler that cleares the resources, when the positions change
        mesh->positionsChanged += [this](vislab::Mesh* mesh, const void*)
        {
            this->clearMeshResources(mesh);
        };

        return resource;
    }

    unsigned int Renderer::getTextureResources(vislab::ColormapTexture* texture)
    {
        auto it = mTextures.find(texture);
        if (it != mTextures.end())
        {
            return it->second;
        }

        if (!texture->scalarField)
            throw std::logic_error("Scalar field needed needed!");

        unsigned int tex;
        glGenTextures(1, &tex);
        glBindTexture(GL_TEXTURE_2D, tex);
        // set the texture wrapping/filtering options (on the currently bound texture object)
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        // load and generate the texture
        Eigen::Vector2i resolution = texture->scalarField->getGrid()->getResolution();

        static std::vector<Eigen::Vector4f> data;
        data.resize(resolution.prod() * 4);
        for (int i = 0; i < resolution.prod(); ++i)
        {
            float value           = texture->scalarField->getArray()->getValue(i).x();
            Eigen::Vector4d color = texture->transferFunction.map(value);
            data[i]       = color.cwiseMin(Eigen::Vector4d::Ones()).cwiseMax(Eigen::Vector4d::Zero()).cast<float>();
        }
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, resolution.x(), resolution.y(), 0, GL_RGBA, GL_FLOAT, data.data());

        // memorize the tuple
        mTextures.insert(std::make_pair(texture, tex));

        // hook an event handler that cleares the resources, when the positions change
        texture->scalarFieldChanged += [this](vislab::ColormapTexture* texture, const void*)
        {
            this->clearTextureResources(texture);
        };

        return tex;
    }
}
