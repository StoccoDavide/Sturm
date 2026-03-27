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

template <typename T>
class SequenceTest : public ::testing::Test {};

using TestTypes = ::testing::Types<float, double>;
TYPED_TEST_SUITE(SequenceTest, TestTypes);

TYPED_TEST(SequenceTest, Test1) {
  using T = TypeParam;

  Poly<T> p1(3);
  p1 << 1.0, -3.0, 2.0;  // p1(x) = 1 - 3x + 2x^2

  Sequence<T> seq(p1);

  typename Poly<T>::Vector sol_seq0(3);
  sol_seq0 << 1.0 / 3.0, -1.0, 2.0 / 3.0;

  EXPECT_TRUE(seq.get(0).coeffs().isApprox(sol_seq0));
  EXPECT_EQ(seq.get(0).degree(), 2);
  EXPECT_EQ(seq.get(0).order(), 3);

  typename Poly<T>::Vector sol_seq1(2);
  sol_seq1 << -3.0 / 4.0, 1.0;

  EXPECT_TRUE(seq.get(1).coeffs().isApprox(sol_seq1));
  EXPECT_EQ(seq.get(1).degree(), 1);
  EXPECT_EQ(seq.get(1).order(), 2);

  typename Poly<T>::Vector sol_seq2(1);
  sol_seq2 << 1.0;

  EXPECT_TRUE(seq.get(2).coeffs().isApprox(sol_seq2));
  EXPECT_EQ(seq.get(2).degree(), 0);
  EXPECT_EQ(seq.get(2).order(), 1);
}

TYPED_TEST(SequenceTest, Test2) {
  using T = TypeParam;

  Poly<T> p2(3);
  p2 << 0.0, 1.0, 1.0;  // p2(x) = x + x^2

  Sequence<T> seq(p2);

  typename Poly<T>::Vector sol_seq0(3);
  sol_seq0 << 0.0, 1.0, 1.0;

  EXPECT_TRUE(seq.get(0).coeffs().isApprox(sol_seq0));
  EXPECT_EQ(seq.get(0).degree(), 2);
  EXPECT_EQ(seq.get(0).order(), 3);

  typename Poly<T>::Vector sol_seq1(2);
  sol_seq1 << 1.0 / 2.0, 1.0;

  EXPECT_TRUE(seq.get(1).coeffs().isApprox(sol_seq1));
  EXPECT_EQ(seq.get(1).degree(), 1);
  EXPECT_EQ(seq.get(1).order(), 2);

  typename Poly<T>::Vector sol_seq2(1);
  sol_seq2 << 1.0;

  EXPECT_TRUE(seq.get(2).coeffs().isApprox(sol_seq2));
  EXPECT_EQ(seq.get(2).degree(), 0);
  EXPECT_EQ(seq.get(2).order(), 1);
}
