# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                    #
#                                                                             #
# The Sandals project is distributed under the GNU GPLv3.                     #
#                                                                             #
# Davide Stocco                                             Enrico Bertolazzi #
# University of Trento                                   University of Trento #
# e-mail: davide.stocco@unitn.it           e-mail: enrico.bertolazzi@unitn.it #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

file(GLOB_RECURSE
  NONLINEAR_SOLVER_TESTS
  "${CMAKE_CURRENT_SOURCE_DIR}/NonlinearSolver/test_*.cc"
)
add_executable(test_nonlinear_solver ${NONLINEAR_SOLVER_TESTS})
target_link_libraries(test_nonlinear_solver PRIVATE sandals)

file(GLOB_RECURSE
RUNGEKUTTA_TESTS
  "${CMAKE_CURRENT_SOURCE_DIR}/RungeKutta/test_*.cc"
)
add_executable(test_rungekutta ${RUNGEKUTTA_TESTS})
target_link_libraries(test_rungekutta PRIVATE sandals)

# That's all folks!
