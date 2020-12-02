# Entry point for user

HEADERS += \
    demographic_cycle.h \
    demographic_sim.h \
    env_grid_cell.h \
    env_param.h \
    environment.h \
    funder_data.h \
    funders.h \
    funders_success.h \
    grid_view.h \
    grn.h \
    ind_param.h \
    individual.h \
    meta_param.h \
    phenotype.h \
    pop_param.h \
    population.h \
    relaxation.hpp \
    sim_parameters.h \
    sim_view.h \
    simulation.h \
    utilities.h

SOURCES += \
    demographic_cycle.cpp \
    demographic_sim.cpp \
    env_grid_cell.cpp \
    env_param.cpp \
    environment.cpp \
    funder_data.cpp \
    funders.cpp \
    funders_success.cpp \
    grid_view.cpp \
    grn.cpp \
    ind_param.cpp \
    individual.cpp \
    main_logic_only.cpp.cpp \
    meta_param.cpp \
    phenotype.cpp \
    pop_param.cpp \
    population.cpp \
    relaxation.cpp \
    sim_parameters.cpp \
    sim_view.cpp \
    simulation.cpp \
    utilities.cpp


INCLUDEPATH += libs
HEADERS += $$PWD/libs/hrtree/*.hpp
