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

TEMPLATE_TEST_CASE("GCD", "[template]", float, double) {

  SECTION("Test 1") {
    Poly<TestType> p1(3); p1 << 1, -3, 2; // p1(x) = 1 - 3x + 2x^2
    Poly<TestType> p2(3); p2 << 0, 1, 1; // p2(x) = x + x^2
    Poly<TestType> gcd; Sturm::GCD<TestType>(p1, p2, gcd);
    Vector<TestType> sol_gcd(1); sol_gcd << 1; // gcd(x) = 1
    REQUIRE(gcd.coeffs().isApprox(sol_gcd));
    REQUIRE(gcd.degree() == 0);
  }

  SECTION("Test 2") {
    Poly<TestType> p3(4); p3 << 1, 0, -1, 0; // p3(x) = 1 - x^2
    Poly<TestType> p4(3); p4 << 0, 0, 1; // p4(x) = x^2
    Poly<TestType> gcd; Sturm::GCD<TestType>(p3, p4, gcd);
    Vector<TestType> sol_gcd(1); sol_gcd << 1.0; // gcd(x) = 1
    REQUIRE(gcd.coeffs().isApprox(sol_gcd));
    REQUIRE(gcd.degree() == 0);
  }

  SECTION("Test 3") {
    Poly<TestType> p5(3); p5 << 1, -2, 1; // p5(x) = 1 - 2x + x^2
    Poly<TestType> p6(2); p6 << 1, -1; // p6(x) = 1 - x
    Poly<TestType> gcd; Sturm::GCD<TestType>(p5, p6, gcd);
    Vector<TestType> sol_gcd(2); sol_gcd << 1, -1; // gcd(x) = 1 - x
    REQUIRE(gcd.coeffs().isApprox(sol_gcd));
    REQUIRE(gcd.degree() == 1);
  }

}
