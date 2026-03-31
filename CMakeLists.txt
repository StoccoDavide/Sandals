# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                    #
#                                                                             #
# The Sandals project is distributed under the BSD 2-Clause License.          #
#                                                                             #
# Davide Stocco                                             Enrico Bertolazzi #
# University of Trento                                   University of Trento #
# davide.stocco@unitn.it                           enrico.bertolazzi@unitn.it #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

cmake_minimum_required(VERSION 3.14)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set(SANDALS_FULL_VERSION "0.0.0")

project(Sandals
  VERSION ${SANDALS_FULL_VERSION}
  DESCRIPTION "A C++ library for ODEs/DAEs integration"
  HOMEPAGE_URL "https://github.com/StoccoDavide/Sandals"
  LANGUAGES CXX
)

message(STATUS "Sandals: version ${SANDALS_FULL_VERSION}")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

if(NOT CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/lib")
  message(STATUS "Sandals: CMAKE_INSTALL_PREFIX not specified, defaulting to ${CMAKE_CURRENT_SOURCE_DIR}")
endif()

set(OPTIMIST_EXTERNAL_EIGEN OFF)
if(TARGET eigen)
  message(STATUS "Optimist: Found external eigen")
  set(OPTIMIST_EXTERNAL_EIGEN ON)
endif()

# Dependencies
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external)

# Library definition
add_library(Sandals INTERFACE)
add_library(Sandals::Sandals ALIAS Sandals)

target_link_libraries(Sandals INTERFACE Eigen3::Eigen Optimist::Optimist)

target_include_directories(Sandals INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/core/include>
  $<INSTALL_INTERFACE:include>
)

# Source files
file(GLOB_RECURSE RELATIVE_HEADER_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  "${CMAKE_CURRENT_SOURCE_DIR}/core/include/*.hh")

# Tests
option(SANDALS_BUILD_TESTS "Build tests" OFF)
if(SANDALS_BUILD_TESTS)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/core/tests)
endif()

# Benchmarks
option(SANDALS_BUILD_BENCHMARKS "Build benchmarks" OFF)
if(SANDALS_BUILD_BENCHMARKS)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/core/benchmarks)
endif()

# Examples
option(SANDALS_BUILD_EXAMPLES "Build examples" OFF)
if(SANDALS_BUILD_EXAMPLES)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/core/examples)
endif()

# Jacobian checks
option(SANDALS_CHECK_JACOBIANS "Check Jacobians" OFF)
if(SANDALS_CHECK_JACOBIANS)
  target_compile_definitions(Sandals INTERFACE SANDALS_CHECK_JACOBIANS)
endif()

# Plotting
option(SANDALS_ENABLE_PLOTTING "Enable plotting" OFF)
if(SANDALS_ENABLE_PLOTTING)
  target_compile_definitions(Sandals INTERFACE SANDALS_ENABLE_PLOTTING)
endif()

# Set compilation flags
if(MSVC)
  target_compile_options(Sandals INTERFACE /W4 /WX)
else()
  target_compile_options(Sandals INTERFACE -Wall -Wextra -pedantic -Werror)
endif()

# Installation
if (NOT SANDALS_EXTERNAL_EIGEN)
  install(
    TARGETS eigen
    EXPORT OptimistTargets
  )
endif()

install(
  TARGETS Sandals
  EXPORT SandalsTargets
  )

install(
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
  DESTINATION include
)

include(CMakePackageConfigHelpers)

export(EXPORT SandalsTargets
  FILE ${CMAKE_CURRENT_BINARY_DIR}/SandalsTargets.cmake
  NAMESPACE Sandals::
)

install(EXPORT SandalsTargets
  FILE SandalsTargets.cmake
  NAMESPACE Sandals::
  DESTINATION lib/cmake/Sandals
)

write_basic_package_version_file(
  SandalsConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/SandalsConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/SandalsConfig.cmake
  INSTALL_DESTINATION lib/cmake/Sandals
)

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/SandalsConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/SandalsVersion.cmake
  DESTINATION lib/cmake/Sandals
)
