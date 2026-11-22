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

# file(GLOB_RECURSE TABLEAU_TESTS "${CMAKE_CURRENT_SOURCE_DIR}/Tableau/test_*.cc")
# add_executable(test_tableau ${TABLEAU_TESTS})
# target_link_libraries(test_tableau PRIVATE Sandals GTest::gtest)

# file(GLOB_RECURSE SYSTEM_TESTS "${CMAKE_CURRENT_SOURCE_DIR}/System/test_*.cc")
# add_executable(test_system ${SYSTEM_TESTS})
# include_directories(${CMAKE_CURRENT_SOURCE_DIR}/DynamicSystems)
# target_link_libraries(test_system PRIVATE Sandals GTest::gtest ROOT::Core ROOT::Gpad)

# file(GLOB_RECURSE ORDER_TESTS "${CMAKE_CURRENT_SOURCE_DIR}/Order/test_*.cc")
# add_executable(test_order ${ORDER_TESTS})
# include_directories(${CMAKE_CURRENT_SOURCE_DIR}/DynamicSystems)
# target_link_libraries(test_order PRIVATE Sandals GTest::gtest)

file(GLOB_RECURSE PROBLEM_TESTS "${CMAKE_CURRENT_SOURCE_DIR}/Problem/test_*.cc")
add_executable(test_problem ${PROBLEM_TESTS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/DynamicSystems)
target_link_libraries(test_problem PRIVATE Sandals GTest::gtest ROOT::Core ROOT::Gpad)

file(GLOB_RECURSE SENSITIVITY_TESTS "${CMAKE_CURRENT_SOURCE_DIR}/Sensitivity/test_*.cc")
add_executable(test_sensitivity ${SENSITIVITY_TESTS})
target_link_libraries(test_sensitivity PRIVATE Sandals GTest::gtest)
add_test(NAME test_sensitivity COMMAND test_sensitivity)

file(GLOB_RECURSE HESSIAN_TESTS "${CMAKE_CURRENT_SOURCE_DIR}/Hessian/test_*.cc")
add_executable(test_hessian ${HESSIAN_TESTS})
target_link_libraries(test_hessian PRIVATE Sandals GTest::gtest)
add_test(NAME test_hessian COMMAND test_hessian)

file(GLOB_RECURSE SHOOTING_TESTS "${CMAKE_CURRENT_SOURCE_DIR}/Shooting/test_*.cc")
add_executable(test_shooting ${SHOOTING_TESTS})
target_link_libraries(test_shooting PRIVATE Sandals GTest::gtest)
add_test(NAME test_shooting COMMAND test_shooting)

add_executable(test_collocation ${CMAKE_CURRENT_SOURCE_DIR}/OptimalControl/test_collocation.cc)
target_include_directories(test_collocation PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/OptimalControl)
target_link_libraries(test_collocation PRIVATE Sandals GTest::gtest)
add_test(NAME test_collocation COMMAND test_collocation)

if(TARGET Sandals::Ipopt)
  add_executable(test_ipopt ${CMAKE_CURRENT_SOURCE_DIR}/OptimalControl/test_ipopt.cc)
  target_include_directories(test_ipopt PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/OptimalControl)
  target_link_libraries(test_ipopt PRIVATE Sandals::Ipopt GTest::gtest)
  add_test(NAME test_ipopt COMMAND test_ipopt WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endif()
