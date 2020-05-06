#ifndef GLSL_PROGRAM_HPP_INCLUCED
#define GLSL_PROGRAM_HPP_INCLUDED

#include <string>
#include <initializer_list>
#include "glsl.hpp"
#include "filesystem.hpp"


namespace glsl {

  GLuint shader_from_str(const char* source, GLenum shaderType);
  GLuint shader_from_file(const fs::path& glslFile, GLenum shaderType);
  std::string shader_source_from_file(const fs::path& glslFile);

  GLuint link_program(std::initializer_list<GLuint> shader, bool deleteShader);

}

#endif
