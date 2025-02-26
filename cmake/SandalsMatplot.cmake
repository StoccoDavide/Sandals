if (SANDALS_ENABLE_PLOTTING)
  set(MATPLOT_REQUIRED_VERSION 1.2.0)
  cmake_policy(SET CMP0135 NEW)

  list(APPEND CMAKE_PREFIX_PATH "${SANDALS_THIRD_PARTY_DIR}")
  find_package(
    Matplot++
    ${MATPLOT_REQUIRED_VERSION}
    NO_MODULE
    QUIET
  )

  get_target_property(MATPLOT_INCLUDE_DIRS
    Matplot++::matplot
    INTERFACE_INCLUDE_DIRECTORIES
  )

  message(STATUS "Sandals: Found Matplot++ installed in ${MATPLOT_INCLUDE_DIRS}")

  if(NOT TARGET Matplot++::matplot)
    message(STATUS "Sandals: Did not find Matplot++ ${MATPLOT_REQUIRED_VERSION} installed, please "
      "install it manually via your package manager.")
  endif()
else()
  add_library(Matplot++::matplot INTERFACE IMPORTED)
  message(STATUS "Sandals: Matplot++ plotting disabled")
endif()
