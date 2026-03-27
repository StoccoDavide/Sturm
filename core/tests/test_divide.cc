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
class DivisionTest : public ::testing::Test {};

using TestTypes = ::testing::Types<float, double>;
TYPED_TEST_SUITE(DivisionTest, TestTypes);

TYPED_TEST(DivisionTest, P1DivP2) {
  using T = TypeParam;

  Poly<T> p1(3);
  p1 << 1.0, -3.0, 2.0;
  Poly<T> p2(3);
  p2 << 0.0, 1.0, 1.0;

  Poly<T> q, r;
  Sturm::divide<T>(p1, p2, q, r);

  typename Poly<T>::Vector sol_quot(1);
  sol_quot << 2.0;

  typename Poly<T>::Vector sol_rem(2);
  sol_rem << 1.0, -5.0;

  EXPECT_TRUE(q.coeffs().isApprox(sol_quot));
  EXPECT_TRUE(r.coeffs().isApprox(sol_rem));
  EXPECT_EQ(q.degree(), 0);
  EXPECT_EQ(r.degree(), 1);
}

TYPED_TEST(DivisionTest, P2DivP1) {
  using T = TypeParam;

  Poly<T> p1(3);
  p1 << 1.0, -3.0, 2.0;
  Poly<T> p2(3);
  p2 << 0.0, 1.0, 1.0;

  Poly<T> q, r;
  Sturm::divide<T>(p2, p1, q, r);

  typename Poly<T>::Vector sol_quot(1);
  sol_quot << 1.0 / 2.0;

  typename Poly<T>::Vector sol_rem(2);
  sol_rem << -1.0 / 2.0, 5.0 / 2.0;

  EXPECT_TRUE(q.coeffs().isApprox(sol_quot));
  EXPECT_TRUE(r.coeffs().isApprox(sol_rem));
  EXPECT_EQ(q.degree(), 0);
  EXPECT_EQ(r.degree(), 1);
}
