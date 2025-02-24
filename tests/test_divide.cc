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

TEMPLATE_TEST_CASE("Division", "[template]", float, double) {

  Poly<TestType> p1(3); p1 << 1.0, -3.0, 2.0; // p1(x) = 1 - 3x + 2x^2
  Poly<TestType> p2(3); p2 << 0.0, 1.0, 1.0; // p2(x) = x + x^2

  SECTION("Test 1") {
    Poly<TestType> q, r; Sturm::divide<TestType>(p1, p2, q, r);
    Vector<TestType> sol_quot(1); sol_quot << 2.0; // q(x) = 2
    Vector<TestType> sol_rem(2); sol_rem << 1.0, -5.0; // r(x) = x - 5
    REQUIRE(q.coeffs().isApprox(sol_quot));
    REQUIRE(r.coeffs().isApprox(sol_rem));
    REQUIRE(q.degree() == 0);
    REQUIRE(r.degree() == 1);
  }

  SECTION("Test 2") {
    Poly<TestType> q, r; Sturm::divide<TestType>(p2, p1, q, r);
    Vector<TestType> sol_quot(1); sol_quot << 1.0/2.0; // q(x) = 1/2
    Vector<TestType> sol_rem(2); sol_rem << -1.0/2.0, 5.0/2.0; // r(x) = -1/2 + 5/2x
    REQUIRE(q.coeffs().isApprox(sol_quot));
    REQUIRE(r.coeffs().isApprox(sol_rem));
    REQUIRE(q.degree() == 0);
    REQUIRE(r.degree() == 1);
  }

}
