/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco, Mattia Piazza and Enrico Bertolazzi.                       *
 *                                                                                               *
 * The Sturm project is distributed under the BSD 2-Clause License.                              *
 *                                                                                               *
 * Davide Stocco                          Mattia Piazza                        Enrico Bertolazzi *
 * University of Trento               University of Trento                  University of Trento *
 * davide.stocco@unitn.it            mattia.piazza@unitn.it           enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "Sturm.hh"

namespace Sturm {

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Sequence::build( Poly const & P ) {
    m_intervals.clear();
    Poly DP, M, R;
    P.derivative( DP );
    this->m_sturm.clear();
    this->m_sturm.reserve(P.order());
    this->m_sturm.emplace_back(P);  this->m_sturm.back().adjust_degree();
    this->m_sturm.emplace_back(DP); this->m_sturm.back().adjust_degree();
    Integer ns = 1;
    while ( true ) {
      divide( this->m_sturm[ns-1], this->m_sturm[ns], M, R );
      if ( R.order() <= 0 ) break;
      this->m_sturm.emplace_back(-R);
      ++ns;
    }
    // divide by GCD
    for ( Integer i = 0; i < ns; ++i ) {
      divide( this->m_sturm[i], this->m_sturm.back(), M, R );
      M.normalize();
      this->m_sturm[i] = M;
    }
    this->m_sturm.back().set_scalar(1);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Integer
  Sequence::sign_variations( Real x, bool & on_root ) const {
    Integer sign_var  = 0;
    Integer last_sign = 0;
    Integer npoly     = Integer(this->m_sturm.size());
    Real    v         = this->m_sturm[0].evaluate(x);
    on_root = false;
    if      ( v > 0 ) last_sign = 1;
    else if ( v < 0 ) last_sign = -1;
    else {
      on_root   = true;
      last_sign = 0;
    }
    for ( Integer i = 1; i < npoly; ++i ) {
      v = this->m_sturm[i].evaluate(x);
      if ( v > 0 ) {
        if ( last_sign == -1 ) ++sign_var;
        last_sign = 1;
      } else if ( v < 0 ) {
        if ( last_sign == 1 ) ++sign_var;
        last_sign = -1;
      }
    }
    return sign_var;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


  Integer
  Sequence::separate_roots( Real a_in, Real b_in ) {
    using std::abs;
    using std::max;

    m_intervals.clear();
    m_intervals.reserve( this->m_sturm.size() );

    Interval I0, I1;
    m_a = I0.a = a_in;
    m_b = I0.b = b_in;

    I0.va = sign_variations( I0.a, I0.a_on_root );
    I0.vb = sign_variations( I0.b, I0.b_on_root );

    Integer n_roots = std::abs( I0.va - I0.vb );

    if ( n_roots <= 1 ) {
      if ( n_roots == 1 && !I0.a_on_root && !I0.b_on_root ) {
        m_intervals.push_back(I0);
      }
      if ( I0.a_on_root ) {
        I1.a  = I1.b  = I0.a;
        I1.va = I1.vb = I0.va;
        I1.a_on_root = I1.b_on_root = true;
        m_intervals.push_back(I1);
      }
      if ( I0.b_on_root ) {
        I1.a  = I1.b  = I0.b;
        I1.va = I1.vb = I0.vb;
        I1.a_on_root = I1.b_on_root = true;
        m_intervals.push_back(I1);
      }
      return Integer(m_intervals.size());
    }

    // search intervals
    std::vector<Interval> I_stack; I_stack.clear(); I_stack.reserve(this->m_sturm.size());
    I_stack.push_back(I0);
    while ( I_stack.size() > 0 ) {
      I0 = I_stack.back(); I_stack.pop_back();
      // controllo se una sola radice
      n_roots = std::abs( I0.va - I0.vb );
      if ( n_roots <= 1 ) {
        if ( I0.a_on_root ) {
          I0.b         = I0.a;
          I0.vb        = I0.va;
          I0.b_on_root = true;
          m_intervals.push_back(I0);
        } else if ( I0.b_on_root ) {
          I0.a         = I0.b;
          I0.va        = I0.vb;
          I0.a_on_root = true;
          m_intervals.push_back(I0);
        } else if ( n_roots == 1 ) {
          m_intervals.push_back(I0);
        }
      } else if ( abs(I0.b-I0.a) <= 10*EPSILON*max(Real(1),max(abs(I0.b),abs(I0.a))) ) {
        I1.a  = I1.b  = I0.a;
        I1.va = I1.vb = 0;
        I1.a_on_root = I1.b_on_root = true;
        I_stack.push_back(I1);
      } else {
        Real    c  = (I0.a+I0.b)/2;
        bool    c_on_root;
        Integer vc = sign_variations( c, c_on_root );
        // check interval [a,c]
        if ( I0.va != vc || c_on_root || I0.a_on_root ) {
          if ( c < I1.b ) { // check if it is a true reduction
            I1.a = I0.a; I1.va = I0.va; I1.a_on_root = I0.a_on_root;
            I1.b = c;    I1.vb = vc;    I1.b_on_root = c_on_root;
            I_stack.push_back(I1);
          } else if ( c_on_root ) {
            I1.a         = I1.b         = c;
            I1.a_on_root = I1.b_on_root = true;
            I1.va        = I1.vb        = 0;
            I_stack.push_back(I1);
          } else if ( I0.a_on_root ) {
            I1.a         = I1.b         = I0.a;
            I1.a_on_root = I1.b_on_root = true;
            I1.va        = I1.vb        = 0;
            I_stack.push_back(I1);
          }
        }
        // check interval [c,b]
        if ( I0.vb != vc || I0.b_on_root ) {
          if ( c > I0.a ) {
            I1.a = c;    I1.va = vc;    I1.a_on_root = c_on_root;
            I1.b = I0.b; I1.vb = I0.vb; I1.b_on_root = I0.b_on_root;
            I_stack.push_back(I1);
          } else if ( I0.b_on_root ) {
            I1.a         = I1.b         = I0.b;
            I1.a_on_root = I1.b_on_root = true;
            I1.va        = I1.vb        = 0;
            I_stack.push_back(I1);
          }
        }
      }
    }
    // sort intervals
    std::sort(
      m_intervals.begin(),
      m_intervals.end(),
      []( Interval const & Sa, Interval const & Sb ) { return Sa.a < Sb.a; }
    );
    return Integer(m_intervals.size());
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Integer
  Sequence::separate_roots() {
    // Cauchy's bounds for roots
    Real an  = this->m_sturm[0].leading_coeff();
    Real bnd = 1+this->m_sturm[0].cwiseAbs().maxCoeff()/abs(an);
    return separate_roots( -bnd, bnd );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Sequence::refine_roots() {
    //FIXMEm_fun.setup( &this->m_sturm[0] );
    m_roots.resize( m_intervals.size() );
    Integer n = 0;
    for ( auto & I : m_intervals ) {
      Real & r = m_roots.coeffRef(n++);
      if      ( I.a_on_root ) r = I.a;
      else if ( I.b_on_root ) r = I.b;
      else {
        //FIXMEr = m_solver.eval( I.a, I.b, &m_fun );
        //FIXMEif ( !m_solver.converged() )
          //FIXMEfmt::print( "Warning: Sequence::refine_roots failed at interval N.{}\n", n );
      }
    }
  }

} // namespace Sturm
