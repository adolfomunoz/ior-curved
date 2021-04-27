project(ivp)
cmake_minimum_required(VERSION 3.10)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})

include(Paths)
include(External)
include(Compiler)


###########################################################################################
# TARGETS
###########################################################################################
add_executable(atmosphere main/atmosphere.cc)
add_dependencies(atmosphere eigen svg-cpp-plot ivp mj2)
target_link_libraries(atmosphere ${svg_cpp_plot_libs})
target_compile_definitions(atmosphere PRIVATE ${svg_cpp_plot_defs})

add_executable(rho_graph main/rho_graph.cc)
add_dependencies(rho_graph eigen svg-cpp-plot ivp mj2)
target_link_libraries(rho_graph ${svg_cpp_plot_libs})
target_compile_definitions(rho_graph PRIVATE ${svg_cpp_plot_defs})
