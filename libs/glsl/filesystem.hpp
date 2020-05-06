#ifndef GLSL_FILESYSTEM_HPP_INCLUDED
#define GLSL_FILESYSTEM_HPP_INCLUDED


#if defined _HAS_CXX17
  #include <filesystem>
#else
  #include <experimental/filesystem>
#endif


namespace glsl {

#if defined _HAS_CXX17
  #if (_MSC_VER >= 1920) || (__GNUC__ > 8)
    namespace fs = std::filesystem;
  #else
    namespace fs = std::experimental::filesystem;
  #endif
#else
  namespace fs = std::experimental::filesystem;
#endif

}

#endif
