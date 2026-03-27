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

#include <limits>

using namespace Sturm;

template <typename T>
class PolyTest : public ::testing::Test {};

using TestTypes = ::testing::Types<float, double>;
TYPED_TEST_SUITE(PolyTest, TestTypes);

TYPED_TEST(PolyTest, Addition) {
  using T = TypeParam;

  Poly<T> p1(3);
  p1 << 1.0, -3.0, 2.0;
  Poly<T> p2(3);
  p2 << 0.0, 1.0, 1.0;

  Poly<T> sum(p1 + p2);

  typename Poly<T>::Vector sol_sum(3);
  sol_sum << 1.0, -2.0, 3.0;

  EXPECT_TRUE(sum.coeffs().isApprox(sol_sum));
  EXPECT_EQ(sum.degree(), 2);
  EXPECT_EQ(sum.order(), 3);
}

TYPED_TEST(PolyTest, Subtraction) {
  using T = TypeParam;

  Poly<T> p1(3);
  p1 << 1.0, -3.0, 2.0;
  Poly<T> p2(3);
  p2 << 0.0, 1.0, 1.0;

  Poly<T> sub(p1 - p2);

  typename Poly<T>::Vector sol_sub(3);
  sol_sub << 1.0, -4.0, 1.0;

  EXPECT_TRUE(sub.coeffs().isApprox(sol_sub));
  EXPECT_EQ(sub.degree(), 2);
  EXPECT_EQ(sub.order(), 3);
}

TYPED_TEST(PolyTest, Multiplication) {
  using T = TypeParam;

  Poly<T> p1(3);
  p1 << 1.0, -3.0, 2.0;
  Poly<T> p2(3);
  p2 << 0.0, 1.0, 1.0;

  Poly<T> mul(p1 * p2);

  typename Poly<T>::Vector sol_mul(5);
  sol_mul << 0.0, 1.0, -2.0, -1.0, 2.0;

  EXPECT_TRUE(mul.coeffs().isApprox(sol_mul));
  EXPECT_EQ(mul.degree(), 4);
  EXPECT_EQ(mul.order(), 5);
}

// ---------------------- Evaluation ----------------------
TYPED_TEST(PolyTest, Evaluation) {
  using T = TypeParam;

  Poly<T> p1(3);
  p1 << 1.0, -3.0, 2.0;
  Poly<T> p2(3);
  p2 << 0.0, 1.0, 1.0;

  T eval1 = p1.evaluate(0.0);
  EXPECT_LT(std::abs(eval1 - T(1.0)), std::numeric_limits<T>::epsilon());

  T eval2 = p2.evaluate(2.0);
  EXPECT_LT(std::abs(eval2 - T(6.0)), std::numeric_limits<T>::epsilon());
}

// ---------------------- Differentiation ----------------------
TYPED_TEST(PolyTest, Differentiation) {
  using T = TypeParam;

  Poly<T> p1(3);
  p1 << 1.0, -3.0, 2.0;

  Poly<T> dif;
  p1.derivative(dif);

  typename Poly<T>::Vector sol_dif(2);
  sol_dif << -3.0, 4.0;

  EXPECT_TRUE(dif.coeffs().isApprox(sol_dif));
  EXPECT_EQ(dif.degree(), 1);
  EXPECT_EQ(dif.order(), 2);
}

// ---------------------- Integration ----------------------
TYPED_TEST(PolyTest, Integration) {
  using T = TypeParam;

  Poly<T> p1(3);
  p1 << 1.0, -3.0, 2.0;

  Poly<T> integ;
  p1.integral(integ);

  typename Poly<T>::Vector sol_integ(4);
  sol_integ << 0.0, 1.0, -3.0 / 2.0, 2.0 / 3.0;

  EXPECT_TRUE(integ.coeffs().isApprox(sol_integ));
  EXPECT_EQ(integ.degree(), 3);
  EXPECT_EQ(integ.order(), 4);
}

// ---------------------- Integration (given constant) ----------------------
TYPED_TEST(PolyTest, IntegrationWithConstant) {
  using T = TypeParam;

  Poly<T> p1(3);
  p1 << 1.0, -3.0, 2.0;

  T c{1.0};
  Poly<T> integ;
  p1.integral(integ, c);

  typename Poly<T>::Vector sol_integ(4);
  sol_integ << c, 1.0, -3.0 / 2.0, 2.0 / 3.0;

  EXPECT_TRUE(integ.coeffs().isApprox(sol_integ));
  EXPECT_EQ(integ.degree(), 3);
  EXPECT_EQ(integ.order(), 4);
}
