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
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/catch_template_test_macros.hpp>

using namespace Sturm;

TEMPLATE_TEST_CASE("Root", "[template]", float, double) {

  SECTION("Test 1") {
    Poly<TestType> p1(3); p1 << 1.0, -3.0, 2.0; // p1(x) = 1 - 3x + 2x^2
    Sequence<TestType> seq(p1); Integer n_roots;
    n_roots = seq.separate_roots(-2.0, -1.0);
    REQUIRE(n_roots == 0);
    n_roots = seq.separate_roots(-1.0, 1.0 + Poly<TestType>::EPSILON);
    REQUIRE(n_roots == 2);
    n_roots = seq.separate_roots(-10.0, 10.0);
    REQUIRE(n_roots == 2);
  }

  SECTION("Test 2") {
    Poly<TestType> p2(3); p2 << -1.0, 1.0, 1.0; // p2(x) = -1 + x + x^2
    Sequence<TestType> seq(p2); Integer n_roots;
    n_roots = seq.separate_roots(-2.0, -1.0);
    REQUIRE(n_roots == 1);
    n_roots = seq.separate_roots(-1.0, 1.0);
    REQUIRE(n_roots == 1);
    n_roots = seq.separate_roots(-10.0, 10.0);
    REQUIRE(n_roots == 2);
  }

}
