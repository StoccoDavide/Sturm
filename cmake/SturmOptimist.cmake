set(OPTIMIST_REQUIRED_VERSION 0.0.0)
cmake_policy(SET CMP0135 NEW)

list(APPEND CMAKE_PATH_PREFIX "${STURM_THIRD_PARTY_DIR}")
find_package(
  Optimist
  ${OPTIMIST_REQUIRED_VERSION}
  NO_MODULE
  QUIET
)

if(NOT TARGET Optimist::Optimist)
  message(STATUS "Sturm: Did not find Optimist ${OPTIMIST_REQUIRED_VERSION} installed, downloading to "
    "${STURM_THIRD_PARTY_DIR}")
  include(FetchContent)

  set(FETCHCONTENT_BASE_DIR "${STURM_THIRD_PARTY_DIR}")
  fetchcontent_declare(
    Optimist
    # URL "https://github.com/StoccoDavide/Optimist/archive/refs/tags/${OPTIMIST_REQUIRED_VERSION}.tar.gz"
    GIT_REPOSITORY "https://github.com/StoccoDavide/Optimist"
    GIT_TAG        main
  )

  fetchcontent_makeavailable(Optimist)

  if(NOT TARGET Optimist::Optimist)
    message(FATAL_ERROR "Optimist::Optimist target was not found after fetchcontent")
  endif()
else()
  get_target_property(OPTIMIST_INCLUDE_DIRS
    Optimist::Optimist
    INTERFACE_INCLUDE_DIRECTORIES
  )
  message(STATUS "Sturm: Found Optimist installed in ${OPTIMIST_INCLUDE_DIRS}")
endif()
