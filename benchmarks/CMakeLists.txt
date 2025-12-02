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

file(GLOB_RECURSE TWPBVPC_RUN "${CMAKE_CURRENT_SOURCE_DIR}/TWPBVPC/run_twpbvpc.cc")
add_executable(run_twpbvpc ${TWPBVPC_RUN})
target_link_libraries(run_twpbvpc PRIVATE Sandals GTest::gtest)

file(GLOB_RECURSE TWPBVPC_PLOT "${CMAKE_CURRENT_SOURCE_DIR}/TWPBVPC/plot_twpbvpc.cc")
add_executable(plot_twpbvpc ${TWPBVPC_PLOT})
target_link_libraries(plot_twpbvpc PRIVATE Sandals ROOT::Core ROOT::Gpad)
