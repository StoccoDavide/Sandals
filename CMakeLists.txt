# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                    #
#                                                                             #
# The Sandals project is distributed under the GNU GPLv3.                     #
#                                                                             #
# Davide Stocco                                             Enrico Bertolazzi #
# University of Trento                                   University of Trento #
# e-mail: davide.stocco@unitn.it           e-mail: enrico.bertolazzi@unitn.it #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

cmake_minimum_required(VERSION 3.10)

project(sandals
  VERSION 0.0.0
  DESCRIPTION "A header-only library to harness the power of BSplines"
  HOMEPAGE_URL "https://github.com/StoccoDavide/Sandals"
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)
set(THIRD_PARTY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/third_party")

option(BUILD_TESTS "Build tests" OFF)
option(BUILD_BENCHMARKS "Build benchmarks" OFF)
option(BUILD_EXAMPLES "Build examples" OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug")
  message("CMAKE_BUILD_TYPE not specified, defaulting to ${CMAKE_BUILD_TYPE}")
endif()

if(NOT CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/lib")
  message("CMAKE_INSTALL_PREFIX not specified, defaulting to ${CMAKE_CURRENT_SOURCE_DIR}")
endif()

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
include_directories("${EIGEN3_INCLUDE_DIR}")

add_library(sandals INTERFACE)

target_include_directories(sandals INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}/include
)
message(STATUS "CMAKE_CURRENT_SOURCE_DIR = ${CMAKE_CURRENT_SOURCE_DIR}")
file(GLOB_RECURSE HEADER_FILES_HH "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hh")
message(STATUS "HEADERS_HH  = ${HEADER_FILES_HH}")
file(GLOB_RECURSE HEADER_FILES_HXX "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hxx")
message(STATUS "HEADERS_HXX = ${HEADER_FILES_HXX}")
target_sources(sandals INTERFACE ${HEADER_FILES_HH} ${HEADER_FILES_HXX})

if(BUILD_TESTS)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/tests)
endif()

if(BUILD_BENCHMARKS)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/benchmarks)
endif()

if(BUILD_EXAMPLES)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/examples)
endif()

if(MSVC)
  target_compile_options(sandals INTERFACE /W4 /WX)
else()
  target_compile_options(sandals INTERFACE -Wall -Wextra -pedantic -Werror)
endif()


# That's all folks!
