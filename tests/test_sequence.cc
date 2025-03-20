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

TEMPLATE_TEST_CASE("Sequence", "[template]", float, double) {

  SECTION("Test 1") {
    Poly<TestType> p1(3); p1 << 1.0, -3.0, 2.0; // p1(x) = 1 - 3x + 2x^2
    Sequence<TestType> seq(p1);
    typename Poly<TestType>::Vector sol_seq0(3); sol_seq0 << 1.0/3.0, -1.0, 2.0/3.0;
    REQUIRE(seq.get(0).coeffs().isApprox(sol_seq0));
    REQUIRE(seq.get(0).degree() == 2);
    REQUIRE(seq.get(0).order() == 3);
    typename Poly<TestType>::Vector sol_seq1(2); sol_seq1 << -3.0/4.0, 1.0;
    REQUIRE(seq.get(1).coeffs().isApprox(sol_seq1));
    REQUIRE(seq.get(1).degree() == 1);
    REQUIRE(seq.get(1).order() == 2);
    typename Poly<TestType>::Vector sol_seq2(1); sol_seq2 << 1.0;
    REQUIRE(seq.get(2).coeffs().isApprox(sol_seq2));
    REQUIRE(seq.get(2).degree() == 0);
    REQUIRE(seq.get(2).order() == 1);
  }

  SECTION("Test 2") {
    Poly<TestType> p2(3); p2 << 0.0, 1.0, 1.0; // p2(x) = x + x^2
    Sequence<TestType> seq(p2);
    typename Poly<TestType>::Vector sol_seq0(3); sol_seq0 << 0.0, 1.0, 1.0;
    REQUIRE(seq.get(0).coeffs().isApprox(sol_seq0));
    REQUIRE(seq.get(0).degree() == 2);
    REQUIRE(seq.get(0).order() == 3);
    typename Poly<TestType>::Vector sol_seq1(2); sol_seq1 << 1.0/2.0, 1.0;
    REQUIRE(seq.get(1).coeffs().isApprox(sol_seq1));
    REQUIRE(seq.get(1).degree() == 1);
    REQUIRE(seq.get(1).order() == 2);
    typename Poly<TestType>::Vector sol_seq2(1); sol_seq2 << 1.0;
    REQUIRE(seq.get(2).coeffs().isApprox(sol_seq2));
    REQUIRE(seq.get(2).degree() == 0);
    REQUIRE(seq.get(2).order() == 1);
  }

}
