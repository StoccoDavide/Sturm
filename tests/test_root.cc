/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sturm project is distributed under the BSD 2-Clause License.                              *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * davide.stocco@unitn.it                                             enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Sturm library
#include "Sturm.hh"

// Catch2 library
#include <catch2/catch_test_macros.hpp>

using namespace Sturm;

TEST_CASE("Root") {

  Poly p1(3); p1 << 1.0, -3.0, 2.0; // p1(x) = 1 - 3x + 2x^2
  Poly p2(3); p2 << 0.0, 1.0, 1.0; // p2(x) = x + x^2

  SECTION("Number") {
    Sequence seq(p1); Integer n_roots;
    n_roots = seq.separate_roots(-2.0, -1.0);
    REQUIRE(n_roots == 0);
    n_roots = seq.separate_roots(-1.0, 1.0);
    REQUIRE(n_roots == 1);
    n_roots = seq.separate_roots(-10.0, 10.0);
    REQUIRE(n_roots == 2);
  }

}
