/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sturm project is distributed under the BSD 2-Clause License.                              *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * davide.stocco@unitn.it                                             enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef INCLUDE_STURM_SEQUENCE_HXX
#define INCLUDE_STURM_SEQUENCE_HXX

namespace Sturm
{

  /**
  * \brief Sturm sequence class.
  *
  * This class implements the Sturm sequence for a given polynomial \f$ p(x) \f$. Such sequence is
  * a sequence of polynomials \f$ p_0(x), p_1(x), \ldots, p_n(x) \f$ and allows to compute the number
  * roots in a given interval \f$ [a, b] \f$.
  */
  class Sequence {
  public:
    using Interval = struct Interval {
      Real    a; /**< Lower bound of the interval. */
      Real    b; /**< Upper bound of the interval. */
      Integer va; /**< Sign of the polynomial at the lower bound. */
      Integer vb; /**< Sign of the polynomial at the upper bound. */
      bool    a_on_root; /**< True if the lower bound is a root. */
      bool    b_on_root; /**< True if the upper bound is a root. */
    }; /**< Interval structure. */

  private:
    using Solver = typename Optimist::ScalarRootFinder::Algo748; /**< Solver class of the Sturm sequence. */
    using Function = Solver::FunctionWrapper; /**< Function class of the Sturm sequence. */

    Solver                m_solver;    /**< Solver class of the Sturm sequence. */
    Function              m_function;  /**< Function class of the Sturm sequence. */
    std::vector<Poly>     m_sequence;  /**< Sturm sequence. */
    std::vector<Interval> m_intervals; /**< Computed intervals. */
    Vector                m_roots;     /**< Computed roots. */
    Real                  m_a{0.0};    /**< Lower bound of the interval containing the roots. */
    Real                  m_b{0.0};    /**< Upper bound of the interval containing the roots. */

  public:
    /**
    * Class constructor for the Sturm sequence.
    */
    Sequence() {}

    /**
    * Class constructor for the Sturm sequence given the polynomial \f$ p(x) \f$.
    * \param[in] p Polynomial.
    */
    Sequence(Poly const & p) {this->build(p);}

    /**
    * Get the lower bound of the interval containing the roots.
    * \return The lower bound of the interval containing the roots.
    */
    Real a() const {return this->m_a;}

    /**
    * Get the upper bound of the interval containing the roots.
    * \return The upper bound of the interval containing the roots.
    */
    Real b() const {return this->m_b;}

    /**
    * Given the polynomial \f$ p(x) \f$ build its Sturm sequence
    */
    void build(Poly const & p);

    /**
    * Get the length of the stored Sturm sequence.
    * \return Length of the stored Sturm sequence.
    */
    Integer length() const {return static_cast<Integer>(this->m_sequence.size());}

    /**
    * Get the \f$ i \f$-th polynomial of the stored Sturm sequence.
    * \return The \f$ i \f$-th polynomial of the stored Sturm sequence.
    */
    Poly const & get(Integer i) const {return this->m_sequence[i];}

    /**
    * Compute the sign variations of the stored Sturm sequence at \f$ x \f$.
    * \return The number of sign variations.
    */
    Integer sign_variations(Real x, bool & on_root) const;

    /**
    * Given an interval \f$ [a, b] \f$ compute the subintervals containing a single root.
    * \return The numbers of intervals (roots) found.
    */
    Integer separate_roots(Real a, Real b);

    /**
    * Compute an interval \f$ [a, b] \f$ that contains all the real roots and compute the subintervals
    * containing a single root.
    * \return The numbers of intervals (roots) found.
    */
    Integer separate_roots();

    /**
    * Get the number of roots found.
    * \return The number of roots found.
    */
    Integer roots_number() const {return static_cast<Integer>(this->m_intervals.size());}

    /**
    * Get the \f$ i \f$-th interval containing a single root.
    * \return The \f$ i \f$-th interval containing a single root.
    */
    Interval const & interval(Integer i) const {return this->m_intervals[i];}

    /**
    * Compute the roots in the intervals after the separation.
    */
    void refine_roots();

    /**
    * Get a vector with the computed roots.
    * \return A vector with the computed roots.
    */
    Vector const & roots() const {return this->m_roots;}

  };

  /**
  * Print the Sturm sequence on an output stream.
  * \param[in] os Output stream.
  * \param[in] s Sturm sequence.
  * \return The output stream.
  */
  std::ostream & operator<< (std::ostream & os, Sequence const & s);

} // namespace Sturm

#endif // INCLUDE_STURM_SEQUENCE_HXX