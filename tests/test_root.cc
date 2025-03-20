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
#include "Sturm/Sequence.hh"

// Catch2 library
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/catch_template_test_macros.hpp>

using namespace Sturm;

// Secant lambda function for the root solver with the signature `bool(Real a, Real b, std::function<Real(Real)> f, Real & x)`
template <typename Real>
auto Secant = [](Real a, Real b, std::function<Real(Real)> f, Real & x) {
  Real tolerance{1.0e-12};
  Real fa{f(a)}, fb{f(b)};
  if (std::abs(fa) < tolerance) {x = a; return true;}
  if (std::abs(fb) < tolerance) {x = b; return true;}
  if (fa*fb > 0) {return false;}
  for (int i{0}; i < 1000; ++i) {
    x = a - fa*(b - a)/(fb - fa);
    Real fx{f(x)};
    if (std::abs(fx) < tolerance) {return true;}
    if (fa*fx < 0) {b = x; fb = fx;}
    else {a = x; fa = fx;}
  }
  return false;
};

TEMPLATE_TEST_CASE("Root", "[template]", float, double) {

  SECTION("Test 1") {
    Poly<TestType> p1(3); p1 << 1.0, -3.0, 2.0; // p1(x) = 1 - 3x + 2x^2
    Sequence<TestType> seq(p1); int n_roots;
    n_roots = seq.separate_roots(-2.0, -1.0);
    REQUIRE(n_roots == 0);
    n_roots = seq.separate_roots(-1.0, 1.0 + Poly<TestType>::EPSILON);
    REQUIRE(n_roots == 2);
    n_roots = seq.separate_roots(-10.0, 10.0);
    REQUIRE(n_roots == 2);
    typename Poly<TestType>::Vector roots = seq.refine_roots(Secant<TestType>);
    for (int i{0}; i < roots.size(); ++i) {
      REQUIRE(std::abs(p1.evaluate(roots(i))) < static_cast<TestType>(1.0e-5));
    }
  }

  SECTION("Test 2") {
    Poly<TestType> p2(3); p2 << -1.0, 1.0, 1.0; // p2(x) = -1 + x + x^2
    Sequence<TestType> seq(p2); int n_roots;
    n_roots = seq.separate_roots(-2.0, -1.0);
    REQUIRE(n_roots == 1);
    n_roots = seq.separate_roots(-1.0, 1.0);
    REQUIRE(n_roots == 1);
    n_roots = seq.separate_roots(-10.0, 10.0);
    REQUIRE(n_roots == 2);
    typename Poly<TestType>::Vector roots = seq.refine_roots(Secant<TestType>);
    REQUIRE(roots.size() == 2);
    for (int i{0}; i < roots.size(); ++i) {
      REQUIRE(std::abs(p2.evaluate(roots(i))) < static_cast<TestType>(1.0e-5));
    }
  }

}
