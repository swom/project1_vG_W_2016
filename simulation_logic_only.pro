# Entry point for user

include(simulation.pri)
SOURCES += \
        main_logic_only.cpp

DEFINES += LOGIC_ONLY

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#pragma GCC diagnostic pop
# QMAKE_CXX += -Wunknown-pragmas


CONFIG += c++17
QMAKE_CXXFLAGS += -std=c++17
CONFIG += resources_big

# High warning levels
# SFML goes bad with -Weffc++
QMAKE_CXXFLAGS += -Wall -Wextra -Wshadow -Wnon-virtual-dtor -pedantic

#Allow to compile Hanno's code
#QMAKE_CXXFLAGS += -fopenmp
#LIBS += -fopenmp

# Hanno's code does not compile cleanly under g++
# A warning is an error
# QMAKE_CXXFLAGS += -Werror

# Debug and release settings
CONFIG += debug_and_release
CONFIG(release, debug|release) {
  DEFINES += NDEBUG

  linux{
    # gprof
    QMAKE_CXXFLAGS += -pg
    QMAKE_LFLAGS += -pg
  }
}


# Qt5
QT += core gui

# SFML, default compiling
# GNU/Linux
unix:!macx {
  #LIBS += -lsfml-graphics -lsfml-window -lsfml-system -lsfml-audio
  DEFINES += ON_LINUX

  CONFIG(debug, debug|release) {
    # gcov
    QMAKE_CXXFLAGS += -fprofile-arcs -ftest-coverage
    LIBS += -lgcov
  }
}


win32{
  #LIBS += -lopenal32              #Dependency
  #LIBS += -lfreetype              #Dependency
  LIBS += -lopengl32              #Dependency
  LIBS += -lgdi32                 #Dependency
  LIBS += -lwinmm                 #Dependency
}




