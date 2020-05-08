# Entry point for user

HEADERS += \
    env_grid_cell.h \
    env_param.h \
    environment.h \
    grid_view.h \
    grn.h \
    ind_param.h \
    individual.h \
    phenotype.h \
    pop_param.h \
    population.h \
    relaxation.hpp \
    sim_parameters.h \
    sim_view.h \
    simulation.h

SOURCES += \
    env_grid_cell.cpp \
    env_param.cpp \
    environment.cpp \
    grid_view.cpp \
    grn.cpp \
    ind_param.cpp \
    individual.cpp \
    main.cpp \
    phenotype.cpp \
    pop_param.cpp \
    population.cpp \
    relaxation.cpp \
    sim_parameters.cpp \
    sim_view.cpp \
    simulation.cpp


INCLUDEPATH += libs
HEADERS += $$PWD/libs/hrtree/*.hpp


CONFIG += c++17
QMAKE_CXXFLAGS += -std=c++17
CONFIG += resources_big

# High warning levels
# SFML goes bad with -Weffc++
QMAKE_CXXFLAGS += -Wall -Wextra -Wshadow -Wnon-virtual-dtor -pedantic

# A warning is an error
#QMAKE_CXXFLAGS += -Werror

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
  LIBS += -lsfml-graphics -lsfml-window -lsfml-system -lsfml-audio

  CONFIG(debug, debug|release) {
    # gcov
    QMAKE_CXXFLAGS += -fprofile-arcs -ftest-coverage
    LIBS += -lgcov
  }
}

win32{
  INCLUDEPATH += C:/Qt/sfml/include
  INCLUDEPATH += D:/Qt/sfml/include
  LIBS += -LC:/Qt/sfml/lib
  LIBS += -LD:/Qt/sfml/lib
  CONFIG(debug, debug|release) {
    LIBS += -lsfml-audio -lsfml-graphics -lsfml-window -lsfml-system
  }
  CONFIG(release, debug|release) {
    LIBS += -lsfml-audio-d -lsfml-graphics-d -lsfml-window-d -lsfml-system-d
  }
  #LIBS += -lopenal32              #Dependency
  #LIBS += -lfreetype              #Dependency
  LIBS += -lopengl32              #Dependency
  LIBS += -lgdi32                 #Dependency
  LIBS += -lwinmm                 #Dependency
}




