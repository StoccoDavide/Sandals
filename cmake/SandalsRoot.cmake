if (SANDALS_ENABLE_PLOTTING)
  set(ROOT_REQUIRED_VERSION 6.34.08)
  cmake_policy(SET CMP0135 NEW)

  # list(APPEND CMAKE_PREFIX_PATH "${SANDALS_THIRD_PARTY_DIR}")
  set(ROOT_REQUIRED_COMPONENTS Core Gpad)
  find_package(
    ROOT
    ${ROOT_REQUIRED_VERSION}
    REQUIRED COMPONENTS ${ROOT_REQUIRED_COMPONENTS}
    NO_MODULE
    QUIET
  )

  get_target_property(ROOT_INCLUDE_DIRS
    ROOT::Core
    INTERFACE_INCLUDE_DIRECTORIES
  )

  message(STATUS "Sandals: Found Root installed in ${ROOT_INCLUDE_DIRS}")

  if(NOT TARGET ROOT::Core OR NOT TARGET ROOT::Gpad)
    message(STATUS "Sandals: Did not find Root ${ROOT_REQUIRED_VERSION} installed, please "
      "install it manually via your package manager.")
  endif()
else()
  add_library(ROOT::Core INTERFACE IMPORTED)
  add_library(ROOT::Gpad INTERFACE IMPORTED)
  message(STATUS "Sandals: Root plotting disabled")
endif()
