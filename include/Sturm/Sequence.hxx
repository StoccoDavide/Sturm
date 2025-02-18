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
    };

    // FIXMEclass Algo748_fun : public Algo748_base_fun<Real> {
    // FIXME  Poly<Real> const * P = nullptr;
    // FIXMEpublic:
    // FIXME  void setup( Poly<Real> const * Pin ) { P = Pin; }
    // FIXME  Real eval( Real x ) const override { return P->eval(x); }
    // FIXME};

  private:

    //FIXME: Algo748<Real>    m_solver;
    //FIXME: Algo748_fun      m_fun;

    std::vector<Poly>   m_sturm;
    std::vector<Interval> m_intervals;
    Vector           m_roots;
    Real             m_a{0};
    Real             m_b{0};

  public:

  Sequence() {}

    //!
    //! Given the polynomial \f$ p(x) \f$ build its Sturm sequence
    //!
    void build( Poly const & p );

    //!
    //! Return the length of the stored Sturm sequence.
    //!
    Integer length() const { return Integer(m_sturm.size()); }

    //!
    //! Return the i-th polynomial of the stored Sturm sequence.
    //!
    Poly const & get( Integer i ) const { return m_sturm[i]; }

    //!
    //! Conpute the sign variations of the stored Sturm sequence at \f$ x \f$.
    //!
    Integer sign_variations( Real x, bool & on_root ) const;

    //!
    //! Given an interval \f$ [a,b] \f$
    //! compute the subintervals containing a single root.
    //! Return the numbers of intervals (roots) found.
    //!
    Integer separate_roots( Real a, Real b );

    //!
    //! Compute an interval \f$ [a,b] \f$ that contains all
    //! the real roots and compute the subintervals containing a single root.
    //! Return the numbers of intervals (roots) found.
    //!
    Integer separate_roots();

    Integer n_roots() const { return Integer(m_intervals.size()); }
    Real a() const { return m_a; }
    Real b() const { return m_b; }
    Interval const & get_interval( Integer i ) const { return m_intervals[i]; }

    //!
    //! compute the roots in the intervals after the separation.
    //!
    void refine_roots();

    //!
    //! return a vector with the computed roots
    //!
    Vector const & roots() const { return m_roots; }

  };

  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  */
  template <typename Char>
  inline
  std::basic_ostream<Char> &
  operator << ( std::basic_ostream<Char> & output, Sequence const & S ) {
    output << "Sturm sequence\n";
    for ( Integer i = 0; i < S.length(); ++i )
      output << "P_" << i << "(x) = " << S.get(i) << '\n';

    Integer n = S.n_roots();
    if ( n > 0 ) {
      //FIXMEfmt::print( output, "roots separation for interval [{},{}]\n", S.a(), S.b() );
      for ( Integer i = 0; i < n; ++i ) {
        //FIXMEtypename Sequence::Interval const & I = S.get_interval( i );
        //FIXMEfmt::print( output, "I = [{}, {}], V = [{}, {}]\n", I.a, I.b, I.va, I.vb );
      }
    }
    return output;
  }


} // namespace Sturm

#endif // INCLUDE_STURM_SEQUENCE_HXX