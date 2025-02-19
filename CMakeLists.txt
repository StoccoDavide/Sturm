# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                        #
#                                                                                                 #
# The Sturm project is distributed under the BSD 2-Clause License.                                #
#                                                                                                 #
# Davide Stocco                                                                 Enrico Bertolazzi #
# University of Trento                                                       University of Trento #
# davide.stocco@unitn.it                                               enrico.bertolazzi@unitn.it #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

cmake_minimum_required(VERSION 3.14)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(SturmVersion)
get_version_from_git()

project(Sturm
  VERSION ${STURM_VERSION}
  DESCRIPTION "Computing Sturm sequences in C++"
  HOMEPAGE_URL "https://github.com/StoccoDavide/Sturm"
  LANGUAGES CXX
)

if(NOT STURM_VERSION_PRERELEASE)
  set(STURM_FULL_VERSION "${STURM_VERSION}")
  write_version_file()
else()
  set(STURM_FULL_VERSION "${STURM_VERSION}-${STURM_VERSION_PRERELEASE}")
endif()
message(STATUS "Sturm: version ${STURM_FULL_VERSION}")

set(STURM_THIRD_PARTY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/third_party")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

if(NOT CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/lib")
  message(STATUS  "Sturm: CMAKE_INSTALL_PREFIX not specified, defaulting to ${CMAKE_CURRENT_SOURCE_DIR}")
endif()

# Dependencies
include(SturmEigen3)
include(SturmOptimist)

# Collect source files
file(GLOB_RECURSE HEADER_FILES "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hh" "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hxx")
file(GLOB_RECURSE SOURCE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cxx")
foreach(HEADER_FILE IN LISTS HEADER_FILES)
file(
  RELATIVE_PATH
  RELATIVE_PATH_HEADER_FILE
  ${CMAKE_CURRENT_SOURCE_DIR} ${HEADER_FILE}
)
list(APPEND RELATIVE_HEADER_FILES ${RELATIVE_PATH_HEADER_FILE})
endforeach()

# Library definition
add_library(Sturm SHARED ${SOURCE_FILES} ${HEADER_FILES})
add_library(Sturm::Sturm ALIAS Sturm)

# Link dependencies
target_link_libraries(Sturm PUBLIC Eigen3::Eigen Optimist::Optimist)

target_include_directories(Sturm PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

# Build options
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug")
    message(STATUS "Sturm: CMAKE_BUILD_TYPE not specified, defaulting to ${CMAKE_BUILD_TYPE}")
  endif()

  option(BUILD_TESTS "Build tests" OFF)
  option(BUILD_BENCHMARKS "Build benchmarks" OFF)
  option(BUILD_EXAMPLES "Build examples" OFF)

  if(BUILD_TESTS)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/tests)
  endif()

  if(BUILD_BENCHMARKS)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/benchmarks)
  endif()

  if(BUILD_EXAMPLES)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/examples)
  endif()
endif()

# Compiler options
if(MSVC)
  target_compile_options(Sturm PRIVATE /W4 /WX)
else()
  target_compile_options(Sturm PRIVATE -Wall -Wextra -pedantic -Werror)
endif()

# Installation
install(TARGETS Sturm EXPORT SturmTargets)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION include)

include(CMakePackageConfigHelpers)

export(EXPORT SturmTargets
  FILE ${CMAKE_CURRENT_BINARY_DIR}/SturmTargets.cmake
  NAMESPACE Sturm::
)

install(EXPORT SturmTargets
  FILE SturmTargets.cmake
  NAMESPACE Sturm::
  DESTINATION lib/cmake/Sturm
)

configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/SturmConfig.cmake
  INSTALL_DESTINATION lib/cmake/Sturm
)

write_basic_package_version_file(
  SturmConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/SturmConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/SturmConfigVersion.cmake
  DESTINATION lib/cmake/Sturm
)