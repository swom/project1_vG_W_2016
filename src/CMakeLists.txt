add_executable(simulation 
               demographic_cycle.cpp 
               demographic_sim.cpp 
               environment.cpp 
               env_grid_cell.cpp 
               env_param.cpp 
               funders.cpp 
               funders_success.cpp 
               funder_data.cpp 
               grn.cpp 
               individual.cpp 
               ind_param.cpp 
               main.cpp
               meta_param.cpp
               phenotype.cpp
               population.cpp
               pop_param.cpp
               relaxation.cpp
               simulation.cpp
               sim_parameters.cpp
               utilities.cpp)
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(simulation PUBLIC ${NLOPT_LIB} OpenMP::OpenMP_CXX ${CMAKE_DL_LIBS})
else()
    target_link_libraries(simulation PUBLIC ${NLOPT_LIB} ${CMAKE_DL_LIBS})
endif()
include_directories(${CMAKE_SOURCE_DIR}/libs/ $(CMAKE_SYSTEM_PREFIX_PATH}/include/))
install(TARGETS simulation DESTINATION bin)
install(TARGETS simulation CONFIGURATIONS Release DESTINATION bin)
