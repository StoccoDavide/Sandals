# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                    #
#                                                                             #
# The Sandals project is distributed under the BSD 2-Clause License.          #
#                                                                             #
# Davide Stocco                                             Enrico Bertolazzi #
# University of Trento                                   University of Trento #
# davide.stocco@unitn.it                           enrico.bertolazzi@unitn.it #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# file(GLOB_RECURSE ALP_RIDER_OCP_RUN "${CMAKE_CURRENT_SOURCE_DIR}/AlpRiderOCP/alp_rider_ocp.cc")
# add_executable(alp_rider_ocp ${ALP_RIDER_OCP_RUN})
# target_link_libraries(alp_rider_ocp PRIVATE Sandals GTest::gtest ROOT::Core ROOT::Gpad)
# target_compile_options(alp_rider_ocp PRIVATE -Wno-error=cpp)

# file(GLOB_RECURSE BASIC_OCP_RUN "${CMAKE_CURRENT_SOURCE_DIR}/BasicOCP/basic_ocp.cc")
# add_executable(basic_ocp ${BASIC_OCP_RUN})
# target_link_libraries(basic_ocp PRIVATE Sandals GTest::gtest ROOT::Core ROOT::Gpad)
# target_compile_options(basic_ocp PRIVATE -Wno-error=cpp)

# file(GLOB_RECURSE BANG_BANG_OCP_RUN "${CMAKE_CURRENT_SOURCE_DIR}/BangBangOCP/bang_bang_ocp.cc")
# add_executable(bang_bang_ocp ${BANG_BANG_OCP_RUN})
# target_link_libraries(bang_bang_ocp PRIVATE Sandals GTest::gtest ROOT::Core ROOT::Gpad)
# target_compile_options(bang_bang_ocp PRIVATE -Wno-error=cpp)

file(GLOB_RECURSE DAE_INDEX_2_OCP_RUN "${CMAKE_CURRENT_SOURCE_DIR}/DAEIndex2OCP/dae_index_2_ocp.cc")
add_executable(dae_index_2_ocp ${DAE_INDEX_2_OCP_RUN})
target_link_libraries(dae_index_2_ocp PRIVATE Sandals GTest::gtest ROOT::Core ROOT::Gpad)
target_compile_options(dae_index_2_ocp PRIVATE -Wno-error=cpp)

# file(GLOB_RECURSE MAX_TRAVEL_OCP_RUN "${CMAKE_CURRENT_SOURCE_DIR}/MaxTravelOCP/max_travel_ocp.cc")
# add_executable(max_travel_ocp ${MAX_TRAVEL_OCP_RUN})
# target_link_libraries(max_travel_ocp PRIVATE Sandals GTest::gtest ROOT::Core ROOT::Gpad)
# target_compile_options(max_travel_ocp PRIVATE -Wno-error=cpp)

file(GLOB_RECURSE PENDULUM_OCP_RUN "${CMAKE_CURRENT_SOURCE_DIR}/PendulumOCP/pendulum_ocp.cc")
add_executable(pendulum_ocp ${PENDULUM_OCP_RUN})
target_link_libraries(pendulum_ocp PRIVATE Sandals GTest::gtest ROOT::Core ROOT::Gpad)
target_compile_options(pendulum_ocp PRIVATE -Wno-error=cpp)

file(GLOB_RECURSE PENDULUM_GGL_OCP_RUN "${CMAKE_CURRENT_SOURCE_DIR}/PendulumGGLOCP/pendulum_ggl_ocp.cc")
add_executable(pendulum_ggl_ocp ${PENDULUM_GGL_OCP_RUN})
target_link_libraries(pendulum_ggl_ocp PRIVATE Sandals GTest::gtest ROOT::Core ROOT::Gpad)
target_compile_options(pendulum_ggl_ocp PRIVATE -Wno-error=cpp)

file(GLOB_RECURSE RAYLEIGH_OCP_RUN "${CMAKE_CURRENT_SOURCE_DIR}/RayleighOCP/rayleigh_ocp.cc")
add_executable(rayleigh_ocp ${RAYLEIGH_OCP_RUN})
target_link_libraries(rayleigh_ocp PRIVATE Sandals GTest::gtest ROOT::Core ROOT::Gpad)
target_compile_options(rayleigh_ocp PRIVATE -Wno-error=cpp)

# file(GLOB_RECURSE TWPBVPC_RUN "${CMAKE_CURRENT_SOURCE_DIR}/TWPBVPC/run_twpbvpc.cc")
# add_executable(run_twpbvpc ${TWPBVPC_RUN})
# target_link_libraries(run_twpbvpc PRIVATE Sandals GTest::gtest)

# file(GLOB_RECURSE TWPBVPC_PLOT "${CMAKE_CURRENT_SOURCE_DIR}/TWPBVPC/plot_twpbvpc.cc")
# add_executable(plot_twpbvpc ${TWPBVPC_PLOT})
# target_link_libraries(plot_twpbvpc PRIVATE Sandals ROOT::Core ROOT::Gpad)
# target_compile_options(plot_twpbvpc PRIVATE -Wno-error=cpp)
