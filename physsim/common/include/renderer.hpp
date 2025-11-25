#pragma once

#include <vislab/opengl/constant_buffer.hpp>

#include <Eigen/Eigen>

#include <memory>
#include <unordered_map>

namespace vislab
{
    class Scene;
    class OpenGL;
    class ShaderOpenGL;
    class Mesh;
    class ColormapTexture;
}

namespace physsim
{
    /**
     * @brief Basic OpenGL renderer for the physsim course.
     */
    class Renderer
    {
    public:
        /**
         * @brief Creates the device resources.
         * @param opengl Reference to the openGL handle.
         * @return True, if the creation succeeded.
         */
        bool createDevice(vislab::OpenGL* opengl);

        /**
         * @brief Creates the swap chain resources, i.e., resources that depend on the screen resolution.
         * @param opengl Reference to the openGL handle.
         * @return True, if the creation succeeded.
         */
        bool createSwapChain(vislab::OpenGL* opengl);

        /**
         * @brief Releases the device resources.
         */
        void releaseDevice();

        /**
         * @brief Releases the swap chain resources.
         */
        void releaseSwapChain();

        /**
         * @brief Draws the frame by invoking OpenGL calls.
         * @param opengl Reference to the OpenGL context.
         */
        void draw(vislab::OpenGL* opengl);

        /**
         * @brief Reference of the scene, which contains all the geometry to render.
         */
        std::shared_ptr<vislab::Scene> scene;

        /**
         * @brief Clears resources that are stored for a certain mesh.
         * @param mesh Mesh to clear resources from.
        */
        void clearMeshResources(vislab::Mesh* mesh);

        /**
         * @brief Clears resources that are stored for a certain texture.
         * @param texture Texture to clear resources from.
         */
        void clearTextureResources(vislab::ColormapTexture* texture);

    protected:
        /**
         * @brief Uniform buffer parameters.
         */
        struct Params
        {
            /**
             * @brief Model world matrix.
             */
            Eigen::Matrix4f worldMatrix;

            /**
             * @brief Matrix for transforming normals to world space.
             */
            Eigen::Matrix4f normalMatrix;

            /**
             * @brief View matrix.
             */
            Eigen::Matrix4f viewMatrix;

            /**
             * @brief Projection matrix.
             */
            Eigen::Matrix4f projMatrix;

            /**
             * @brief Inverse of model view projection matrix.
             */
            Eigen::Matrix4f invMVPMatrix;

            /**
             * @brief World position of the camera.
             */
            Eigen::Vector4f eyePosition;

            /**
             * @brief Diffuse color.
             */
            Eigen::Vector4f diffuse;

            /**
             * @brief Point light positions.
             */
            Eigen::Vector4f lightPos[8];

            /**
             * @brief Point light intensities.
             */
            Eigen::Vector4f lightIntensity[8];
        };

        /**
         * @brief Wrapper that holds a uniform buffer.
         */
        vislab::ConstantBufferOpenGL<Params> mParams;

        /**
         * @brief Shader program for the rendering of spheres.
         */
        std::shared_ptr<vislab::ShaderOpenGL> mSphereShader;

        /**
         * @brief Shader program for the rendering of rectangles.
         */
        std::shared_ptr<vislab::ShaderOpenGL> mRectangleShader;

        /**
         * @brief Shader program for the rendering of rectangles with texture.
         */
        std::shared_ptr<vislab::ShaderOpenGL> mTexturedRectangleShader;

        /**
         * @brief Shader program for the rendering of triangles with flat shading.
         */
        std::shared_ptr<vislab::ShaderOpenGL> mFlatTriangleShader;

        /**
         * @brief Shader program for the rendering of triangles with phong shading.
         */
        std::shared_ptr<vislab::ShaderOpenGL> mTriangleShader;

        /**
         * @brief Vertex array object for the sphere.
         */
        uint32_t mSphereVAO;

        /**
         * @brief Vertex array object for the rectangle.
         */
        uint32_t mRectangleVAO;

        /**
         * @brief Collection of OpenGL handles related to a triangle mesh.
        */
        struct MeshResource
        {
            /**
             * @brief Vertex buffer of positions.
            */
            unsigned int vbo_positions;

            /**
             * @brief Vertex buffer of normals.
            */
            unsigned int vbo_normals;

            /**
             * @brief Index buffer.
            */
            unsigned int ibo;

            /**
             * @brief Vertex array object.
            */
            unsigned int vao;
        };

        /**
         * @brief Helper function that does lazy construction of VBOs, IBOs, and VAOs from a mesh.
         * @param mesh Mesh to map to GL resources.
         * @return Handle to the mesh resources.
        */
        MeshResource getMeshResources(vislab::Mesh* mesh);

        /**
         * @brief Collection of VBOs IBOs, and VAOs for the triangle meshes.
        */
        std::unordered_map<vislab::Mesh*, MeshResource> mTriangleMeshes;

        /**
         * @brief Helper function that does lazy construction of textures.
         * @param texture Color mapped texture to build GL resources for.
         * @return Handle to texture.
         */
        unsigned int getTextureResources(vislab::ColormapTexture* texture);

        /**
         * @brief Collection of texture handles for color mapped textures.
         */
        std::unordered_map<vislab::ColormapTexture*, unsigned int> mTextures;
    };
}
