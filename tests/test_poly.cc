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

TEMPLATE_TEST_CASE("Poly", "[template]", float, double) {

  Poly<TestType> p1(3); p1 << 1.0, -3.0, 2.0; // p1(x) = 1 - 3x + 2x^2
  Poly<TestType> p2(3); p2 << 0.0, 1.0, 1.0; // p2(x) = x + x^2

  SECTION("Addition") {
    Poly<TestType> sum(p1 + p2); Vector<TestType> sol_sum(3); sol_sum << 1.0, -2.0, 3.0;
    REQUIRE(sum.coeffs().isApprox(sol_sum)); // sum(x) = 1 - 2x + 3x^2
    REQUIRE(sum.degree() == 2);
    REQUIRE(sum.order() == 3);
  }

  SECTION("Subtraction") {
    Poly<TestType> sub(p1 - p2); Vector<TestType> sol_sub(3); sol_sub << 1.0, -4.0, 1.0;
    REQUIRE(sub.coeffs().isApprox(sol_sub)); // sub(x) = 1 - 4x + x^2
    REQUIRE(sub.degree() == 2);
    REQUIRE(sub.order() == 3);
  }

  SECTION("Multiplication") {
    Poly<TestType> mul(p1 * p2); Vector<TestType> sol_mul(5); sol_mul << 0.0, 1.0, -2.0, -1.0, 2.0;
    REQUIRE(mul.coeffs().isApprox(sol_mul)); // mul(x) = x + 2x^2 - x^3 - 2x^4
    REQUIRE(mul.degree() == 4);
    REQUIRE(mul.order() == 5);
  }

  SECTION("Evaluation") {
    double eval1{p1.evaluate(0.0)};
    REQUIRE((eval1 - double(1.0)) < std::numeric_limits<TestType>::epsilon()); // p1(1) = 1
    double eval2{p2.evaluate(2.0)};
    REQUIRE((eval2 - double(6.0)) < std::numeric_limits<TestType>::epsilon()); // p2(2) = 6
  }

  SECTION("Differentiation") {
    Poly<TestType> dif; p1.derivative(dif);
    Vector<TestType> sol_dif(2); sol_dif << -3.0, 4.0;
    REQUIRE(dif.coeffs().isApprox(sol_dif)); // dif(x) = -3 + 4x
    REQUIRE(dif.degree() == 1);
    REQUIRE(dif.order() == 2);
  }

  SECTION("Integration") {
    Poly<TestType> integ; p1.integral(integ);
    Vector<TestType> sol_integ(4); sol_integ << 0.0, 1.0, -3.0/2.0, 2.0/3.0;
    REQUIRE(integ.coeffs().isApprox(sol_integ)); // integ(x) = x - 3/2x^2 + 2/3x^3
    REQUIRE(integ.degree() == 3);
    REQUIRE(integ.order() == 4);
  }

  SECTION("Integration (given constant)") {
    double c{1.0}; Poly<TestType> integ; p1.integral(integ, c);
    Vector<TestType> sol_integ(4); sol_integ << c, 1.0, -3.0/2.0, 2.0/3.0;
    REQUIRE(integ.coeffs().isApprox(sol_integ)); // integ(x) = c + x - 3/2x^2 + 2/3x^3
    REQUIRE(integ.degree() == 3);
    REQUIRE(integ.order() == 4);
  }

}
