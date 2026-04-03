# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                    #
#                                                                             #
# The Sandals project is distributed under the BSD 2-Clause License.          #
#                                                                             #
# Davide Stocco                                             Enrico Bertolazzi #
# University of Trento                                   University of Trento #
# davide.stocco@unitn.it                           enrico.bertolazzi@unitn.it #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

file(GLOB_RECURSE CURVE_OCP_RUN "${CMAKE_CURRENT_SOURCE_DIR}/CurveOCP/curve_ocp.cc")
add_executable(curve_ocp ${CURVE_OCP_RUN})
target_link_libraries(curve_ocp PRIVATE Sandals GTest::gtest ROOT::Core ROOT::Gpad)
target_compile_options(curve_ocp PRIVATE -Wno-error=cpp)

file(GLOB_RECURSE PENDULUM_OCP_RUN "${CMAKE_CURRENT_SOURCE_DIR}/PendulumOCP/pendulum_ocp.cc")
add_executable(pendulum_ocp ${PENDULUM_OCP_RUN})
target_link_libraries(pendulum_ocp PRIVATE Sandals GTest::gtest ROOT::Core ROOT::Gpad)
target_compile_options(pendulum_ocp PRIVATE -Wno-error=cpp)

file(GLOB_RECURSE TWPBVPC_RUN "${CMAKE_CURRENT_SOURCE_DIR}/TWPBVPC/run_twpbvpc.cc")
add_executable(run_twpbvpc ${TWPBVPC_RUN})
target_link_libraries(run_twpbvpc PRIVATE Sandals GTest::gtest)

file(GLOB_RECURSE TWPBVPC_PLOT "${CMAKE_CURRENT_SOURCE_DIR}/TWPBVPC/plot_twpbvpc.cc")
add_executable(plot_twpbvpc ${TWPBVPC_PLOT})
target_link_libraries(plot_twpbvpc PRIVATE Sandals ROOT::Core ROOT::Gpad)
target_compile_options(plot_twpbvpc PRIVATE -Wno-error=cpp)
