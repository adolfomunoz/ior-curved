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

add_executable(zenith_rho_graph main/zenith_rho_graph.cc)
add_dependencies(zenith_rho_graph eigen svg-cpp-plot ivp mj2)
target_link_libraries(zenith_rho_graph ${svg_cpp_plot_libs})
target_compile_definitions(zenith_rho_graph PRIVATE ${svg_cpp_plot_defs})

add_executable(rho_ch_cone_graph main/rho_ch_cone_graph.cc)
add_dependencies(rho_ch_cone_graph eigen svg-cpp-plot ivp mj2)
target_link_libraries(rho_ch_cone_graph ${svg_cpp_plot_libs})
target_compile_definitions(rho_ch_cone_graph PRIVATE ${svg_cpp_plot_defs})

add_executable(3d_cone_graph main/3d_cone_graph.cc)
add_dependencies(3d_cone_graph eigen svg-cpp-plot ivp mj2)
target_link_libraries(3d_cone_graph ${svg_cpp_plot_libs})
target_compile_definitions(3d_cone_graph PRIVATE ${svg_cpp_plot_defs})

add_executable(atmosphere_cone main/atmosphere_cone.cc)
add_dependencies(atmosphere_cone eigen svg-cpp-plot ivp mj2)
target_link_libraries(atmosphere_cone ${svg_cpp_plot_libs})
target_compile_definitions(atmosphere_cone PRIVATE ${svg_cpp_plot_defs})

add_executable(height_graph main/height_graph.cc)
add_dependencies(height_graph eigen svg-cpp-plot ivp mj2)
target_link_libraries(height_graph ${svg_cpp_plot_libs})
target_compile_definitions(height_graph PRIVATE ${svg_cpp_plot_defs})

add_executable(trajectory main/trajectory.cc)
add_dependencies(trajectory eigen svg-cpp-plot ivp mj2 patterns rapidxml)
target_link_libraries(trajectory ${svg_cpp_plot_libs})
target_compile_definitions(trajectory PRIVATE ${svg_cpp_plot_defs})

add_executable(ior main/ior.cc)
add_dependencies(ior eigen svg-cpp-plot ivp mj2 patterns rapidxml)
target_link_libraries(ior ${svg_cpp_plot_libs})
target_compile_definitions(ior PRIVATE ${svg_cpp_plot_defs})

