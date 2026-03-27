# # # # # # # # # # # # # # # # # # # # # # #  # # # # # # # # # # # # # # # # #
# Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                     #
#                                                                              #
# The Sturm project is distributed under the BSD 2-Clause License.             #
#                                                                              #
# Davide Stocco                                              Enrico Bertolazzi #
# University of Trento                                    University of Trento #
# davide.stocco@unitn.it                            enrico.bertolazzi@unitn.it #
# # # # # # # # # # # # # # # # # # # # # # # #  # # # # # # # # # # # # # # # #

cmake_minimum_required(VERSION 3.14)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

project(Sturm
  VERSION 0.0.0
  DESCRIPTION "Computing Sturm sequences in C++"
  HOMEPAGE_URL "https://github.com/StoccoDavide/Sturm"
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

if(NOT CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/lib")
  message(STATUS "Sturm: CMAKE_INSTALL_PREFIX not specified, defaulting to ${CMAKE_CURRENT_SOURCE_DIR}")
endif()

# Dependencies
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external)

# Library definition
add_library(Sturm INTERFACE)
add_library(Sturm::Sturm ALIAS Sturm)

target_link_libraries(Sturm INTERFACE Eigen3::Eigen)

target_include_directories(Sturm INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/core/include>
  $<INSTALL_INTERFACE:include>
)

# Source files
file(GLOB_RECURSE HEADER_FILES_HH "${CMAKE_CURRENT_SOURCE_DIR}/core/include/*.hh")
foreach(HEADER_FILE IN LISTS HEADER_FILES_HH)
  file(
    RELATIVE_PATH
    RELATIVE_PATH_HEADER_FILE
    ${CMAKE_CURRENT_SOURCE_DIR} ${HEADER_FILE}
  )
  list(APPEND RELATIVE_HEADER_FILES ${RELATIVE_PATH_HEADER_FILE})
endforeach()
target_sources(Sturm INTERFACE $<BUILD_INTERFACE:${RELATIVE_HEADER_FILES}>)

# Tests
option(STURM_BUILD_TESTS "Build tests" OFF)
if(STURM_BUILD_TESTS)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/core/tests)
endif()

# Set compilation flags
if(MSVC)
  target_compile_options(Sturm INTERFACE /W4 /WX)
else()
  target_compile_options(Sturm INTERFACE -Wall -Wextra -pedantic -Werror)
endif()

# Installation
install(
  TARGETS Sturm
  EXPORT SturmTargets
  DESTINATION lib
)

install(
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/core/include/
  DESTINATION include
)

# export(EXPORT SturmTargets
#   FILE ${CMAKE_CURRENT_BINARY_DIR}/SturmTargets.cmake
#   NAMESPACE Sturm::
# )

# install(EXPORT SturmTargets
#   FILE SturmTargets.cmake
#   NAMESPACE Sturm::
#   DESTINATION lib/cmake/Sturm
# )

# configure_package_config_file(
#   ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Config.cmake.in
#   ${CMAKE_CURRENT_BINARY_DIR}/SturmConfig.cmake
#   INSTALL_DESTINATION lib/cmake/Sturm
# )

# write_basic_package_version_file(
#   SturmConfigVersion.cmake
#   VERSION ${PROJECT_VERSION}
#   COMPATIBILITY SameMajorVersion
# )

# install(FILES
#   ${CMAKE_CURRENT_BINARY_DIR}/SturmConfig.cmake
#   ${CMAKE_CURRENT_BINARY_DIR}/SturmConfigVersion.cmake
#   DESTINATION lib/cmake/Sturm
# )
