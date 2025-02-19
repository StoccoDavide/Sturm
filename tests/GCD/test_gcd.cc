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

TEST_CASE("GCD") {

  SECTION("Test 1") {
    Poly p1(3); p1 << 1, -3, 2; // p1(x) = 1 - 3x + 2x^2
    Poly p2(3); p2 << 0, 1, 1; // p2(x) = x + x^2
    Poly gcd; Sturm::GCD(p1, p2, gcd);
    Vector sol_gcd(1); sol_gcd << 1; // gcd(x) = 1
    REQUIRE(gcd.coeffs().isApprox(sol_gcd));
    REQUIRE(gcd.degree() == 0);
  }

  SECTION("Test 2") {
    Poly p3(4); p3 << 1, 0, -1, 0; // p3(x) = 1 - x^2
    Poly p4(3); p4 << 0, 0, 1; // p4(x) = x^2
    Poly gcd; Sturm::GCD(p3, p4, gcd);
    Vector sol_gcd(1); sol_gcd << 1.0; // gcd(x) = 1
    REQUIRE(gcd.coeffs().isApprox(sol_gcd));
    REQUIRE(gcd.degree() == 0);
  }

  SECTION("Test 3") {
    Poly p5(3); p5 << 1, -2, 1; // p5(x) = 1 - 2x + x^2
    Poly p6(2); p6 << 1, -1; // p6(x) = 1 - x
    Poly gcd; Sturm::GCD(p5, p6, gcd);
    Vector sol_gcd(2); sol_gcd << 1, -1; // gcd(x) = 1 - x
    REQUIRE(gcd.coeffs().isApprox(sol_gcd));
    REQUIRE(gcd.degree() == 1);
  }

}
