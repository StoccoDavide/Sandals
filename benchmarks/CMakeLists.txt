# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                        #
#                                                                                                 #
# The Sandals project is distributed under the BSD 2-Clause License.                              #
#                                                                                                 #
# Davide Stocco                                                                 Enrico Bertolazzi #
# University of Trento                                                       University of Trento #
# e-mail: davide.stocco@unitn.it                               e-mail: enrico.bertolazzi@unitn.it #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

include(SandalsGTest)
include(SandalsRoot)

file(GLOB_RECURSE TWPBVPC_PLOTS "${CMAKE_CURRENT_SOURCE_DIR}/TWPBVPC/plot_twpbvpc.cc")
add_executable(plot_twpbvpc ${TWPBVPC_PLOTS})
target_link_libraries(plot_twpbvpc PRIVATE Sandals ROOT::Core ROOT::Gpad)

file(GLOB_RECURSE TWPBVPC_BENCHMARK "${CMAKE_CURRENT_SOURCE_DIR}/TWPBVPC/benchmark_twpbvpc.cc")
add_executable(benchmark_twpbvpc ${TWPBVPC_BENCHMARK})
target_link_libraries(benchmark_twpbvpc PRIVATE Sandals GTest::gtest)
