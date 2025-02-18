# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Copyright (c) 2025, Davide Stocco, Mattia Piazza and Enrico Bertolazzi.                         #
#                                                                                                 #
# The Sturm project is distributed under the BSD 2-Clause License.                                #
#                                                                                                 #
# Davide Stocco                           Mattia Piazza                         Enrico Bertolazzi #
# University of Trento                University of Trento                   University of Trento #
# davide.stocco@unitn.it             mattia.piazza@unitn.it            enrico.bertolazzi@unitn.it #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

include(SturmCatch2)

file(GLOB_RECURSE TEST_POLY "${CMAKE_CURRENT_SOURCE_DIR}/Poly/test_poly.cc")
add_executable(test_poly ${TEST_POLY})
target_link_libraries(test_poly PRIVATE Sturm)

file(GLOB_RECURSE TEST_STURM "${CMAKE_CURRENT_SOURCE_DIR}/Sturm/test_sturm.cc")
add_executable(test_sturm ${TEST_STURM})
target_link_libraries(test_sturm PRIVATE Sturm)

file(GLOB_RECURSE TEST_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/Root/test_root.cc")
add_executable(test_root ${TEST_ROOT})
target_link_libraries(test_root PRIVATE Sturm)
