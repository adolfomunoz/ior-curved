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
add_dependencies(atmosphere svg-cpp-plot ivp mj2 eigen)
target_link_libraries(atmosphere ${svg_cpp_plot_libs})
target_compile_definitions(atmosphere PRIVATE ${svg_cpp_plot_defs})
