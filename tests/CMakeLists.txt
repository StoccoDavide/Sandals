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

file(GLOB_RECURSE
  BROYDEN_TESTS
  "${CMAKE_CURRENT_SOURCE_DIR}/Broyden/test_*.cc"
)
add_executable(test_broyden ${BROYDEN_TESTS})
target_link_libraries(test_broyden PRIVATE sandals Matplot++::matplot)

file(GLOB_RECURSE
  NEWTON_TESTS
  "${CMAKE_CURRENT_SOURCE_DIR}/Newton/test_*.cc"
)
add_executable(test_newton ${NEWTON_TESTS})
target_link_libraries(test_newton PRIVATE sandals Matplot++::matplot)

file(GLOB_RECURSE
RUNGEKUTTA_TESTS
  "${CMAKE_CURRENT_SOURCE_DIR}/RungeKutta/test_*.cc"
)
add_executable(test_rungekutta ${RUNGEKUTTA_TESTS})
target_link_libraries(test_rungekutta PRIVATE sandals Matplot++::matplot)

# That's all folks!
