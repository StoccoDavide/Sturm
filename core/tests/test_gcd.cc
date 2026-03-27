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
#include "Sturm/Poly.hh"

// Google Test
#include <gtest/gtest.h>

using namespace Sturm;
// Sturm library
#include "Sturm/Sequence.hh"

// Google Test
#include <gtest/gtest.h>

using namespace Sturm;

template <typename T>
class GCDTest : public ::testing::Test {};

using TestTypes = ::testing::Types<float, double>;
TYPED_TEST_SUITE(GCDTest, TestTypes);

TYPED_TEST(GCDTest, CoprimePolynomials1) {
  using T = TypeParam;

  Poly<T> p1(3);
  p1 << 1, -3, 2;
  Poly<T> p2(3);
  p2 << 0, 1, 1;

  Poly<T> gcd;
  Sturm::GCD<T>(p1, p2, gcd);

  typename Poly<T>::Vector sol_gcd(1);
  sol_gcd << 1;

  EXPECT_TRUE(gcd.coeffs().isApprox(sol_gcd));
  EXPECT_EQ(gcd.degree(), 0);
}

TYPED_TEST(GCDTest, CoprimePolynomials2) {
  using T = TypeParam;

  Poly<T> p3(4);
  p3 << 1, 0, -1, 0;
  Poly<T> p4(3);
  p4 << 0, 0, 1;

  Poly<T> gcd;
  Sturm::GCD<T>(p3, p4, gcd);

  typename Poly<T>::Vector sol_gcd(1);
  sol_gcd << 1.0;

  EXPECT_TRUE(gcd.coeffs().isApprox(sol_gcd));
  EXPECT_EQ(gcd.degree(), 0);
}

TYPED_TEST(GCDTest, NonTrivialGCD) {
  using T = TypeParam;

  Poly<T> p5(3);
  p5 << 1, -2, 1;
  Poly<T> p6(2);
  p6 << 1, -1;

  Poly<T> gcd;
  Sturm::GCD<T>(p5, p6, gcd);

  typename Poly<T>::Vector sol_gcd(2);
  sol_gcd << 1, -1;

  EXPECT_TRUE(gcd.coeffs().isApprox(sol_gcd));
  EXPECT_EQ(gcd.degree(), 1);
}
