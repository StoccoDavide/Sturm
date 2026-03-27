/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Sturm project is distributed under the BSD 2-Clause License.          *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Sturm library
#include "Sturm/Sequence.hh"

// Google Test
#include <gtest/gtest.h>

using namespace Sturm;

// Secant lambda function for the root solver with the signature `bool(Real a,
// Real b, std::function<Real(Real)> f, Real & x)`
template <typename Real>
auto Secant = [](Real a, Real b, std::function<Real(Real)> f, Real &x) {
  Real tolerance{1.0e-12};
  Real fa{f(a)}, fb{f(b)};
  if (std::abs(fa) < tolerance) {
    x = a;
    return true;
  }
  if (std::abs(fb) < tolerance) {
    x = b;
    return true;
  }
  if (fa * fb > 0) {
    return false;
  }
  for (int i{0}; i < 1000; ++i) {
    x = a - fa * (b - a) / (fb - fa);
    Real fx{f(x)};
    if (std::abs(fx) < tolerance) {
      return true;
    }
    if (fa * fx < 0) {
      b  = x;
      fb = fx;
    } else {
      a  = x;
      fa = fx;
    }
  }
  return false;
};

template <typename T>
class RootTest : public ::testing::Test {};

using TestTypes = ::testing::Types<float, double>;
TYPED_TEST_SUITE(RootTest, TestTypes);

TYPED_TEST(RootTest, Test1) {
  using T = TypeParam;

  Poly<T> p(3);
  p << 1.0, -3.0, 2.0;  // p(x) = 1 - 3x + 2x^2

  Sequence<T> seq(p);
  int n_roots;

  n_roots = seq.separate_roots(-2.0, -1.0);
  EXPECT_EQ(n_roots, 0);

  n_roots = seq.separate_roots(-1.0, 1.0 + Poly<T>::EPSILON);
  EXPECT_EQ(n_roots, 2);

  n_roots = seq.separate_roots(-10.0, 10.0);
  EXPECT_EQ(n_roots, 2);

  typename Poly<T>::Vector roots = seq.refine_roots(Secant<T>);
  EXPECT_EQ(roots.size(), 2);

  for (int i{0}; i < roots.size(); ++i) {
    EXPECT_LT(std::abs(p.evaluate(roots(i))), static_cast<T>(1.0e-5));
  }
}

TYPED_TEST(RootTest, Test2) {
  using T = TypeParam;

  Poly<T> p(3);
  p << -1.0, 1.0, 1.0;  // p(x) = -1 + x + x^2

  Sequence<T> seq(p);
  int n_roots;

  n_roots = seq.separate_roots(-2.0, -1.0);
  EXPECT_EQ(n_roots, 1);

  n_roots = seq.separate_roots(-1.0, 1.0);
  EXPECT_EQ(n_roots, 1);

  n_roots = seq.separate_roots(-10.0, 10.0);
  EXPECT_EQ(n_roots, 2);

  typename Poly<T>::Vector roots = seq.refine_roots(Secant<T>);
  EXPECT_EQ(roots.size(), 2);

  for (int i{0}; i < roots.size(); ++i) {
    EXPECT_LT(std::abs(p.evaluate(roots(i))), static_cast<T>(1.0e-5));
  }
}
