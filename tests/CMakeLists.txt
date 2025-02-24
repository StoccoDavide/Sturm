# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                        #
#                                                                                                 #
# The Sturm project is distributed under the BSD 2-Clause License.                                #
#                                                                                                 #
# Davide Stocco                                                                 Enrico Bertolazzi #
# University of Trento                                                       University of Trento #
# davide.stocco@unitn.it                                               enrico.bertolazzi@unitn.it #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

include(SturmCatch2)

file(GLOB_RECURSE TEST_POLY "${CMAKE_CURRENT_SOURCE_DIR}/test_poly.cc")
add_executable(test_Poly ${TEST_POLY})
target_link_libraries(test_Poly PRIVATE Sturm Catch2::Catch2WithMain)

file(GLOB_RECURSE TEST_DIVIDE "${CMAKE_CURRENT_SOURCE_DIR}/test_divide.cc")
add_executable(test_divide ${TEST_DIVIDE})
target_link_libraries(test_divide PRIVATE Sturm Catch2::Catch2WithMain)

file(GLOB_RECURSE TEST_GCD "${CMAKE_CURRENT_SOURCE_DIR}/test_gcd.cc")
add_executable(test_gcd ${TEST_GCD})
target_link_libraries(test_gcd PRIVATE Sturm Catch2::Catch2WithMain)

file(GLOB_RECURSE TEST_SEQUENCE "${CMAKE_CURRENT_SOURCE_DIR}/test_sequence.cc")
add_executable(test_sequence ${TEST_SEQUENCE})
target_link_libraries(test_sequence PRIVATE Sturm Catch2::Catch2WithMain)

file(GLOB_RECURSE TEST_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/test_root.cc")
add_executable(test_root ${TEST_ROOT})
target_link_libraries(test_root PRIVATE Sturm Catch2::Catch2WithMain)