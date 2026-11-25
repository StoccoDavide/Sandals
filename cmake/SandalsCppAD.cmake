# Optional CppAD (automatic differentiation) support, used by the DAE optimal control problems with
# automatic derivatives (Sandals/OptimalControlAD.hh). Defines the INTERFACE target Sandals::CppAD.
#
# An installed CppAD is found through pkg-config (set PKG_CONFIG_PATH if needed). Otherwise CppAD is
# downloaded, built, and installed once in SANDALS_THIRD_PARTY_DIR at configure time.

option(SANDALS_ENABLE_CPPAD "Enable the optional CppAD automatic differentiation support" ON)

set(SANDALS_CPPAD_TAG "20250000.3")

if(SANDALS_ENABLE_CPPAD)
  find_package(PkgConfig QUIET)

  if(PkgConfig_FOUND)
    pkg_check_modules(CPPAD QUIET IMPORTED_TARGET cppad)
  endif()

  if(NOT CPPAD_FOUND AND PkgConfig_FOUND)
    set(SANDALS_CPPAD_PREFIX "${SANDALS_THIRD_PARTY_DIR}/cppad-install")
    if(NOT EXISTS "${SANDALS_CPPAD_PREFIX}/lib/pkgconfig/cppad.pc")
      message(STATUS "Sandals: CppAD not found, building ${SANDALS_CPPAD_TAG} in ${SANDALS_CPPAD_PREFIX}")
      include(FetchContent)
      set(FETCHCONTENT_BASE_DIR "${SANDALS_THIRD_PARTY_DIR}")
      fetchcontent_declare(
        cppad
        GIT_REPOSITORY "https://github.com/coin-or/CppAD.git"
        GIT_TAG        ${SANDALS_CPPAD_TAG}
        GIT_SHALLOW    TRUE
        SOURCE_SUBDIR  "do-not-add-subdirectory" # download only, CppAD is built below
      )
      fetchcontent_makeavailable(cppad)
      set(SANDALS_CPPAD_BUILD "${SANDALS_THIRD_PARTY_DIR}/cppad-build")
      execute_process(
        COMMAND ${CMAKE_COMMAND} -S ${cppad_SOURCE_DIR} -B ${SANDALS_CPPAD_BUILD}
          -DCMAKE_BUILD_TYPE=Release -Dcppad_prefix=${SANDALS_CPPAD_PREFIX} -Dcppad_cxx_flags=-std=c++20
        RESULT_VARIABLE SANDALS_CPPAD_RESULT OUTPUT_QUIET
      )
      if(SANDALS_CPPAD_RESULT EQUAL 0)
        execute_process(
          COMMAND ${CMAKE_COMMAND} --build ${SANDALS_CPPAD_BUILD} --target install --parallel
          RESULT_VARIABLE SANDALS_CPPAD_RESULT OUTPUT_QUIET
        )
      endif()
      if(NOT SANDALS_CPPAD_RESULT EQUAL 0)
        message(WARNING "Sandals: CppAD build failed, the automatic differentiation support is disabled")
      endif()
    endif()
    set(SANDALS_SAVED_PKG_CONFIG_PATH "$ENV{PKG_CONFIG_PATH}")
    set(ENV{PKG_CONFIG_PATH} "${SANDALS_CPPAD_PREFIX}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
    pkg_check_modules(CPPAD QUIET IMPORTED_TARGET cppad)
    set(ENV{PKG_CONFIG_PATH} "${SANDALS_SAVED_PKG_CONFIG_PATH}")
  endif()

  if(CPPAD_FOUND)
    message(STATUS "Sandals: Found CppAD ${CPPAD_VERSION} in ${CPPAD_INCLUDE_DIRS}")
    add_library(Sandals_CppAD INTERFACE)
    add_library(Sandals::CppAD ALIAS Sandals_CppAD)
    target_link_libraries(Sandals_CppAD INTERFACE Sandals PkgConfig::CPPAD)
  else()
    message(STATUS "Sandals: CppAD not available, the automatic differentiation support is disabled")
  endif()
endif()
