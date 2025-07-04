# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                        #
#                                                                                                 #
# The Sandals project is distributed under the BSD 2-Clause License.                              #
#                                                                                                 #
# Davide Stocco                                                                 Enrico Bertolazzi #
# University of Trento                                                       University of Trento #
# e-mail: davide.stocco@unitn.it                               e-mail: enrico.bertolazzi@unitn.it #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

cmake_minimum_required(VERSION 3.14)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(SandalsVersion)
get_version_from_git()

project(Sandals
  VERSION ${SANDALS_VERSION}
  DESCRIPTION "A C++ library for ODEs/DAEs integration"
  HOMEPAGE_URL "https://github.com/StoccoDavide/Sandals"
  LANGUAGES CXX
)

if(NOT SANDALS_VERSION_PRERELEASE)
  set(SANDALS_FULL_VERSION "${SANDALS_VERSION}")
  write_version_file()
else()
  set(SANDALS_FULL_VERSION "${SANDALS_VERSION}-${SANDALS_VERSION_PRERELEASE}")
endif()
message(STATUS "Sandals: version ${SANDALS_FULL_VERSION}")

set(SANDALS_THIRD_PARTY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/third_party")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

if(NOT CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/lib")
  message(STATUS "Sandals: CMAKE_INSTALL_PREFIX not specified, defaulting to ${CMAKE_CURRENT_SOURCE_DIR}")
endif()

# Dependencies
include(SandalsEigen3)
include(SandalsOptimist)

# Library definition
add_library(Sandals INTERFACE)
add_library(Sandals::Sandals ALIAS Sandals)

target_link_libraries(Sandals INTERFACE Eigen3::Eigen Optimist::Optimist)

target_include_directories(Sandals INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

# INFO: Ensure this options are only available to developers. In theory we could leave them available,
# but the call to `target_sources` breaks `fetchcontent` compatibility due to absolute paths being
# added to `INTERFACE_SOURCES`. I tried solving it, but it seems to be poorly documented, supported, etc.
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug")
    message(STATUS "Sandals: CMAKE_BUILD_TYPE not specified, defaulting to ${CMAKE_BUILD_TYPE}")
  endif()

  option(SANDALS_BUILD_TESTS "Build tests" OFF)
  option(SANDALS_BUILD_BENCHMARKS "Build benchmarks" OFF)
  option(SANDALS_BUILD_EXAMPLES "Build examples" OFF)
  option(SANDALS_ENABLE_PLOTTING "Enable plotting" OFF)

  file(GLOB_RECURSE HEADER_FILES_HH "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hh")
  foreach(HEADER_FILE IN LISTS HEADER_FILES_HH)
    file(
      RELATIVE_PATH
      RELATIVE_PATH_HEADER_FILE
      ${CMAKE_CURRENT_SOURCE_DIR} ${HEADER_FILE}
    )
    list(APPEND RELATIVE_HEADER_FILES ${RELATIVE_PATH_HEADER_FILE})
  endforeach()
  target_sources(Sandals INTERFACE $<BUILD_INTERFACE:${RELATIVE_HEADER_FILES}>)

  if(SANDALS_BUILD_TESTS)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/tests)
  endif()

  if(SANDALS_BUILD_BENCHMARKS)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/benchmarks)
  endif()

  if(SANDALS_BUILD_EXAMPLES)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/examples)
  endif()

  if(SANDALS_ENABLE_PLOTTING)
    target_compile_definitions(Sandals INTERFACE SANDALS_ENABLE_PLOTTING)
  endif()
endif()

if(MSVC)
  target_compile_options(Sandals INTERFACE /W4 /WX)
else()
  target_compile_options(Sandals INTERFACE -Wall -Wextra -pedantic -Werror)
endif()

# Installation
install(TARGETS Sandals EXPORT SandalsTargets)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION include)

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
  ${CMAKE_CURRENT_BINARY_DIR}/SandalsConfigVersion.cmake
  DESTINATION lib/cmake/Sandals
)
