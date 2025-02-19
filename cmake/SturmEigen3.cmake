set(EIGEN_REQUIRED_VERSION 3.4.0)
cmake_policy(SET CMP0135 NEW)

list(APPEND CMAKE_PATH_PREFIX "${STURM_THIRD_PARTY_DIR}")
find_package(
  Eigen3
  ${EIGEN_REQUIRED_VERSION}
  NO_MODULE
  QUIET
)

if(NOT TARGET Eigen3::Eigen)
  message(STATUS "Sturm: Did not find Eigen3 ${EIGEN_REQUIRED_VERSION} installed, downloading to "
    "${STURM_THIRD_PARTY_DIR}")

  include(FetchContent)
  set(FETCHCONTENT_BASE_DIR "${STURM_THIRD_PARTY_DIR}")
  fetchcontent_declare(
    Eigen3
    URL "https://gitlab.com/libeigen/eigen/-/archive/${EIGEN_REQUIRED_VERSION}/eigen-${EIGEN_REQUIRED_VERSION}.tar.gz")

  fetchcontent_makeavailable(Eigen3)

  if(NOT TARGET Eigen3::Eigen)
    message(FATAL_ERROR "Eigen3::Eigen target was not found after fetchcontent")
  endif()
else()
  get_target_property(EIGEN_INCLUDE_DIRS
    Eigen3::Eigen
    INTERFACE_INCLUDE_DIRECTORIES
  )
  message(STATUS "Sturm: Found Eigen3 installed in ${EIGEN_INCLUDE_DIRS}")
endif()
