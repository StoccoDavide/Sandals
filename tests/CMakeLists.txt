# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                    #
#                                                                             #
# The Sandals project is distributed under the GNU GPLv3.                     #
#                                                                             #
# Davide Stocco                                             Enrico Bertolazzi #
# University of Trento                                   University of Trento #
# e-mail: davide.stocco@unitn.it           e-mail: enrico.bertolazzi@unitn.it #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

find_package(Matplot++ REQUIRED)
find_package(GTest REQUIRED)

file(GLOB_RECURSE BROYDEN_TESTS "${CMAKE_CURRENT_SOURCE_DIR}/Broyden/test_*.cc")
add_executable(test_broyden ${BROYDEN_TESTS})
target_link_libraries(test_broyden PRIVATE sandals GTest::gtest)

file(GLOB_RECURSE NEWTON_TESTS "${CMAKE_CURRENT_SOURCE_DIR}/Newton/test_*.cc")
add_executable(test_newton ${NEWTON_TESTS})
target_link_libraries(test_newton PRIVATE sandals GTest::gtest)


file(GLOB_RECURSE EXPLICIT_TESTS "${CMAKE_CURRENT_SOURCE_DIR}/Explicit/test_*.cc")
add_executable(test_explicit ${EXPLICIT_TESTS})
target_link_libraries(test_explicit PRIVATE sandals GTest::gtest Matplot++::matplot)

file(GLOB_RECURSE IMPLICIT_TESTS "${CMAKE_CURRENT_SOURCE_DIR}/Implicit/test_*.cc")
add_executable(test_implicit ${IMPLICIT_TESTS})
target_link_libraries(test_implicit PRIVATE sandals GTest::gtest Matplot++::matplot)

file(GLOB_RECURSE TABLEAU_TESTS "${CMAKE_CURRENT_SOURCE_DIR}/Tableau/test_*.cc")
add_executable(test_tableau ${TABLEAU_TESTS})
target_link_libraries(test_tableau PRIVATE sandals GTest::gtest)

# That's all folks!
