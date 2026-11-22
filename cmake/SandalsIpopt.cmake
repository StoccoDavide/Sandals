# Optional IPOPT interface (header Sandals/Ipopt.hh), found through pkg-config.
# Defines the INTERFACE target Sandals::Ipopt, which links Sandals and IPOPT, when IPOPT is found.

option(SANDALS_ENABLE_IPOPT "Enable the optional IPOPT interface" ON)

if(SANDALS_ENABLE_IPOPT)
  find_package(PkgConfig QUIET)
  if(PkgConfig_FOUND)
    pkg_check_modules(IPOPT QUIET IMPORTED_TARGET ipopt)
    if(NOT IPOPT_FOUND)
      # Fall back to the default MacPorts and Homebrew prefixes
      set(SANDALS_SAVED_PKG_CONFIG_PATH "$ENV{PKG_CONFIG_PATH}")
      set(ENV{PKG_CONFIG_PATH} "$ENV{PKG_CONFIG_PATH}:/opt/local/lib/pkgconfig:/opt/homebrew/lib/pkgconfig:/usr/local/lib/pkgconfig")
      pkg_check_modules(IPOPT QUIET IMPORTED_TARGET ipopt)
      set(ENV{PKG_CONFIG_PATH} "${SANDALS_SAVED_PKG_CONFIG_PATH}")
    endif()
  endif()

  if(IPOPT_FOUND)
    message(STATUS "Sandals: Found IPOPT ${IPOPT_VERSION} in ${IPOPT_INCLUDE_DIRS}")
    add_library(Sandals_Ipopt INTERFACE)
    add_library(Sandals::Ipopt ALIAS Sandals_Ipopt)
    target_link_libraries(Sandals_Ipopt INTERFACE Sandals PkgConfig::IPOPT)
  else()
    message(STATUS "Sandals: IPOPT not found, the IPOPT interface is disabled (install it, e.g., "
      "'sudo port install ipopt' or 'brew install ipopt', or set PKG_CONFIG_PATH)")
  endif()
endif()
