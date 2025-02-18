/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco, Mattia Piazza and Enrico Bertolazzi.                       *
 *                                                                                               *
 * The Sturm project is distributed under the BSD 2-Clause License.                              *
 *                                                                                               *
 * Davide Stocco                          Mattia Piazza                        Enrico Bertolazzi *
 * University of Trento               University of Trento                  University of Trento *
 * davide.stocco@unitn.it            mattia.piazza@unitn.it           enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef INCLUDE_STURM_SEQUENCE_HXX
#define INCLUDE_STURM_SEQUENCE_HXX

namespace Sturm
{


  class Sequence {
  public:

    using Interval = struct Interval {
      Real    a;
      Real    b;
      Integer va;
      Integer vb;
      bool    a_on_root;
      bool    b_on_root;
    }; /**< Interval structure. */

    // FIXMEclass Algo748_fun : public Algo748_base_fun<Real> {
    // FIXME  Poly<Real> const * P = nullptr;
    // FIXMEpublic:
    // FIXME  void setup( Poly<Real> const * Pin ) { P = Pin;}
    // FIXME  Real eval( Real x ) const override {return P->eval(x);}
    // FIXME};

  private:

    //FIXME: Algo748<Real>    m_solver;
    //FIXME: Algo748_fun      m_fun;

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