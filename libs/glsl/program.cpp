#include <vector>
#include <fstream>
#include <iostream>
#include "program.hpp"


namespace glsl {


  namespace {


    void do_collect_includes(const fs::path& glslFile, std::vector<fs::path>& includes)
    {
      auto parent = glslFile.parent_path();
      std::ifstream is(glslFile);
      for (; is;) {
        std::string line;
        std::getline(is, line);
        if (0 == line.compare(0, 10, "#include <")) {
          auto pos = line.rfind('>');
          if (pos != std::string::npos) {
            auto incl = parent / line.substr(10, pos - 10);
            // check for circular includes
            if (incl != glslFile &&
              (includes.cend() == std::find(includes.cbegin(), includes.cend(), incl)))
            {
              do_collect_includes(incl, includes);
            }
          }
        }
      }
      includes.emplace_back(glslFile);
    }


    std::vector<fs::path> collect_includes(const fs::path& glslFile)
    {
      std::vector<fs::path> includes;
      do_collect_includes(glslFile, includes);
      return includes;
    }

  }


  GLuint shader_from_str(const char* source, GLenum shaderType)
  {
    const auto sh = glCreateShader(shaderType);
    glShaderSource(sh, 1, &source, 0);
    glCompileShader(sh);
    return sh;
  }


  std::string shader_source_from_file(const fs::path& glslFile)
  {
    if (!fs::exists(glslFile)) {
      throw std::runtime_error(glslFile.string() + " doesn't exist");
    }
    auto includes = collect_includes(glslFile);
    std::string shader;
    for (auto it = includes.cbegin(); it != includes.cend(); ++it) {
      std::ifstream is(*it);
      for (; is;) {
        std::string line;
        std::getline(is, line);
        if (0 == line.compare(0, 10, "#include <")) {
          shader += "//";   // out-comment include statement
        }
        shader += line;
        shader += '\n';
      }
    }
    return shader;
  }


  GLuint shader_from_file(const fs::path& glslFile, GLenum shaderType)
  {
    auto shader = shader_source_from_file(glslFile);
    return shader_from_str(shader.c_str(), shaderType);
  }


  GLuint link_program(std::initializer_list<GLuint> shader, bool deleteShader)
  {
    GLuint prog = glCreateProgram();
    for (GLuint sh : shader) {
      glAttachShader(prog, sh);
    }
    glLinkProgram(prog);
    for (GLuint sh : shader) {
      glDetachShader(prog, sh);
      if (deleteShader) glDeleteShader(sh);
    }
    GLchar buf[1024];
    GLsizei length;
    glGetProgramInfoLog(prog, 1024, &length, buf);
    if (length) {
#ifdef GLSL_DEBUG_OUTPUT
      glDebugMessageInsert(GL_DEBUG_SOURCE_APPLICATION, GL_DEBUG_TYPE_ERROR, 1, GL_DEBUG_SEVERITY_HIGH, length, buf);
#endif
      glDeleteProgram(prog);
      prog = GL_NONE;
    }
    return prog;
  }

}
