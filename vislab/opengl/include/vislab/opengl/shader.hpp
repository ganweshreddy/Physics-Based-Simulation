#pragma once

#include <string>

namespace vislab
{
    /**
     * @brief Represents an OpenGL shader.
     */
    class ShaderOpenGL
    {
    public:
        /**
         * @brief Constructor.
         */
        ShaderOpenGL();

        /**
         * @brief Destructor.
         */
        ~ShaderOpenGL();

        /**
         * @brief Creates the resources that are independent of the viewport resolution.
         * @return
         */
        bool createDevice();

        /**
         * @brief Releases the resources that are independent of the viewport resolution.
         */
        void releaseDevice();

        /**
         * @brief Bind the shader.
         */
        void bind();

        /**
         * @brief Sets the vertex shader source code.
         * @param source GLSL source code.
         */
        void setVertexShaderSource(const std::string& source);

        /**
         * @brief Sets the geometry shader source code.
         * @param source GLSL source code.
         */
        void setGeometryShaderSource(const std::string& source);

        /**
         * @brief Sets the fragment shader source code.
         * @param source GLSL source code.
         */
        void setPixelShaderSource(const std::string& source);

        /**
         * @brief Gets the OpenGL program handle.
         * @return Program handle.
         */
        unsigned int getProgram() const;

    private:
        /**
         * @brief Compiles a shader.
         * @param shaderSource Shader source to compile.
         * @param shaderType Shader type.
         * @param error Optional error message.
         * @return Shader handle.
         */
        static unsigned int compileShader(const char* shaderSource, unsigned int shaderType, std::string* error = NULL);

        /**
         * @brief Link the shaders into a program.
         * @param error Optional error message.
         * @return Program handle.
         */
        unsigned int linkProgram(std::string* error = NULL);

        /**
         * @brief GLSL vertex shader source code.
         */
        std::string mVertexSource;

        /**
         * @brief GLSL geometry shader source code.
         */
        std::string mGeometrySource;

        /**
         * @brief Fragment/pixel shader source code.
         */
        std::string mPixelSource;

        /**
         * @brief Vertex shader handle.
         */
        unsigned int mVertexShader;

        /**
         * @brief Geometry shader handle.
         */
        unsigned int mGeometryShader;

        /**
         * @brief Fragement/pixel shader handle.
         */
        unsigned int mPixelShader;

        /**
         * @brief Program handle.
         */
        unsigned int mProgram;
    };
}
