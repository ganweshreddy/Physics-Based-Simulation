#include <iostream>
#include <vislab/opengl/opengl.hpp>
#include <vislab/opengl/shader.hpp>

namespace vislab
{
    ShaderOpenGL::ShaderOpenGL()
        : mVertexSource()
        , mGeometrySource()
        , mPixelSource()
        , mVertexShader(0)
        , mGeometryShader(0)
        , mPixelShader(0)
        , mProgram(0)
    {
    }

    ShaderOpenGL::~ShaderOpenGL()
    {
        releaseDevice();
    }

    bool ShaderOpenGL::createDevice()
    {
        std::string error;

        // compile vertex shader
        if (!mVertexSource.empty())
        {
            mVertexShader = compileShader(mVertexSource.c_str(), GL_VERTEX_SHADER, &error);
            if (!error.empty())
                std::cout << "ERROR: Vertex shader compilation failed.\n"
                          << error << std::endl;
        }

        // compile geometry shader
        if (!mGeometrySource.empty())
        {
            mGeometryShader = compileShader(mGeometrySource.c_str(), GL_GEOMETRY_SHADER, &error);
            if (!error.empty())
                std::cout << "ERROR: Geometry shader compilation failed.\n"
                          << error << std::endl;
        }

        // compile pixel shader
        if (!mPixelSource.empty())
        {
            mPixelShader = compileShader(mPixelSource.c_str(), GL_FRAGMENT_SHADER, &error);
            if (!error.empty())
                std::cout << "ERROR: Pixel shader compilation failed.\n"
                          << error << std::endl;
        }

        // link program
        mProgram = linkProgram(&error);
        if (!error.empty())
            std::cout << "ERROR: Shader linkage failed.\n"
                      << error << std::endl;

        // delete the shaders
        if (!mVertexSource.empty())
        {
            glDeleteShader(mVertexShader);
            mVertexShader = 0;
        }
        if (!mGeometrySource.empty())
        {
            glDeleteShader(mGeometryShader);
            mGeometryShader = 0;
        }
        if (!mPixelSource.empty())
        {
            glDeleteShader(mPixelShader);
            mPixelShader = 0;
        }

        return mProgram;
    }

    void ShaderOpenGL::releaseDevice()
    {
        glDeleteProgram(mProgram);
        mProgram = 0;
    }

    void ShaderOpenGL::bind()
    {
        glUseProgram(mProgram);
    }

    void ShaderOpenGL::setVertexShaderSource(const std::string& source)
    {
        mVertexSource = source;
    }

    void ShaderOpenGL::setGeometryShaderSource(const std::string& source)
    {
        mGeometrySource = source;
    }

    void ShaderOpenGL::setPixelShaderSource(const std::string& source)
    {
        mPixelSource = source;
    }

    unsigned int ShaderOpenGL::getProgram() const
    {
        return mProgram;
    }

    unsigned int ShaderOpenGL::compileShader(const char* shaderSource, unsigned int shaderType, std::string* error)
    {
        unsigned int handle = glCreateShader(shaderType);
        glShaderSource(handle, 1, &shaderSource, NULL);
        glCompileShader(handle);
        int success;
        glGetShaderiv(handle, GL_COMPILE_STATUS, &success);
        if (!success && error)
        {
            int infologLength = 0;
            glGetShaderiv(handle, GL_INFO_LOG_LENGTH, &infologLength);
            error->resize(infologLength);
            glGetShaderInfoLog(handle, infologLength, NULL, error->data());
        }
        return handle;
    }

    unsigned int ShaderOpenGL::linkProgram(std::string* error)
    {
        unsigned int handle = glCreateProgram();
        if (mVertexShader)
            glAttachShader(handle, mVertexShader);
        if (mGeometryShader)
            glAttachShader(handle, mGeometryShader);
        if (mPixelShader)
            glAttachShader(handle, mPixelShader);
        glLinkProgram(handle);
        int success;
        glGetProgramiv(handle, GL_LINK_STATUS, &success);
        if (!success && error)
        {
            int infologLength = 0;
            glGetProgramiv(handle, GL_INFO_LOG_LENGTH, &infologLength);
            error->resize(infologLength);
            glGetProgramInfoLog(handle, infologLength, NULL, error->data());
        }
        return handle;
    }
}
