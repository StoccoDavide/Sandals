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
include(SandalsMatplot)

file(GLOB_RECURSE SYSTEM_TESTS "${CMAKE_CURRENT_SOURCE_DIR}/System/test_*.cc")
add_executable(test_system ${SYSTEM_TESTS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/DynamicSystems)
target_link_libraries(test_system PRIVATE Sandals GTest::gtest Matplot++::matplot)

file(GLOB_RECURSE ORDER_TESTS "${CMAKE_CURRENT_SOURCE_DIR}/Order/test_*.cc")
add_executable(test_order ${ORDER_TESTS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/DynamicSystems)
target_link_libraries(test_order PRIVATE Sandals GTest::gtest Matplot++::matplot)

file(GLOB_RECURSE TABLEAU_TESTS "${CMAKE_CURRENT_SOURCE_DIR}/Tableau/test_*.cc")
add_executable(test_tableau ${TABLEAU_TESTS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/DynamicSystems)
target_link_libraries(test_tableau PRIVATE Sandals GTest::gtest)
