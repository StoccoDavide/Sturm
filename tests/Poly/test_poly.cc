/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco, Mattia Piazza and Enrico Bertolazzi.                       *
 *                                                                                               *
 * The Sturm project is distributed under the BSD 2-Clause License.                              *
 *                                                                                               *
 * Davide Stocco                          Mattia Piazza                        Enrico Bertolazzi *
 * University of Trento               University of Trento                  University of Trento *
 * davide.stocco@unitn.it            mattia.piazza@unitn.it           enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Sturm library
#include "Sturm.hh"

// Catch2 library
#include <catch2/catch_test_macros.hpp>

using namespace Sturm;

TEST_CASE("Poly", "[polynomial]") {
  Poly p1(3);
  p1 << 1, -3, 2; // p1(x) = 1 - 3x + 2x^2
  Poly p2(3);
  p2 << 0, 1, 1;  // p2(x) = x + x^2

  SECTION("Addition") {
    Poly sum(p1 + p2);
    Vector sol_sum(3); sol_sum << 1, -2, 3;
    REQUIRE((sum.coeffs() - sol_sum).norm() < EPSILON); // sum(x) = 1 - 2x + 3x^2
  }

  SECTION("Subtraction") {
    Poly diff(p1 - p2);
    Vector sol_diff(3); sol_diff << 1, -4, 1;
    REQUIRE((diff.coeffs() - sol_diff).norm() < EPSILON); // diff(x) = 1 - 4x + x^2
  }

  SECTION("Multiplication") {
    Poly prod(p1 * p2);
    Vector sol_prod(5); sol_prod << 0, 1, -2, -1, 2;
    REQUIRE((prod.coeffs() - sol_prod).norm() < EPSILON); // prod(x) = x + 2x^2 - x^3 - 2x^4
  }

  SECTION("Evaluation") {
    Real eval1{p1.evaluate(0.0)};
    REQUIRE((eval1 - 1.0) < EPSILON); // p1(1) = 1
    Real eval2{p2.evaluate(2.0)};
    REQUIRE((eval2 - 6.0) < EPSILON); // p2(2) = 6
  }

  SECTION("Derivative") {
    Poly deriv;
    p1.derivative(deriv);
    Vector sol_deriv(2); sol_deriv << -3, 4;
    REQUIRE((deriv.coeffs() - sol_deriv).norm() < EPSILON); // deriv(x) = -3 + 4x
  }

  SECTION("Integral") {
    Poly integ;
    p1.integral(integ);
    Vector sol_integ(4); sol_integ << 0, 1, -3.0/2.0, 2.0/3.0;
    REQUIRE((integ.coeffs() - sol_integ).norm() < EPSILON); // integ(x) = x - 3/2x^2 + 2/3x^3
  }
}

