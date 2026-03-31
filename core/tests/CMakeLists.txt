# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                    #
#                                                                             #
# The Sandals project is distributed under the BSD 2-Clause License.          #
#                                                                             #
# Davide Stocco                                             Enrico Bertolazzi #
# University of Trento                                   University of Trento #
# davide.stocco@unitn.it                           enrico.bertolazzi@unitn.it #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

file(GLOB_RECURSE TABLEAU_TESTS "${CMAKE_CURRENT_SOURCE_DIR}/Tableau/test_*.cc")
add_executable(test_tableau ${TABLEAU_TESTS})
target_link_libraries(test_tableau PRIVATE Sandals GTest::gtest)

# file(GLOB_RECURSE SYSTEM_TESTS "${CMAKE_CURRENT_SOURCE_DIR}/System/test_*.cc")
# add_executable(test_system ${SYSTEM_TESTS})
# include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Systems)
# target_link_libraries(test_system PRIVATE Sandals GTest::gtest ROOT::Core ROOT::Gpad)

file(GLOB_RECURSE ORDER_TESTS "${CMAKE_CURRENT_SOURCE_DIR}/Order/test_*.cc")
add_executable(test_order ${ORDER_TESTS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Systems)
target_link_libraries(test_order PRIVATE Sandals GTest::gtest)

# file(GLOB_RECURSE PROBLEM_TESTS "${CMAKE_CURRENT_SOURCE_DIR}/Problem/test_*.cc")
# add_executable(test_problem ${PROBLEM_TESTS})
# include_directories(${CMAKE_CURRENT_SOURCE_DIR}/BoundaryValueProblems)
# target_link_libraries(test_problem PRIVATE Sandals GTest::gtest ROOT::Core ROOT::Gpad)
