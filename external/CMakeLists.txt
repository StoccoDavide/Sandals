# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                    #
#                                                                             #
# The Sandals project is distributed under the BSD 2-Clause License.          #
#                                                                             #
# Davide Stocco                                             Enrico Bertolazzi #
# University of Trento                                   University of Trento #
# davide.stocco@unitn.it                           enrico.bertolazzi@unitn.it #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Eigen
if (NOT SANDALS_EXTERNAL_EIGEN)
  set(EIGEN_BUILD_TESTS OFF)
  add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/eigen)
endif()

# GoogleTest
if (SANDALS_BUILD_TESTS)
  add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/googletest)
endif()
