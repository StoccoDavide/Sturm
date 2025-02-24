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
  * \tparam Real Scalar number type.
  */
  template <typename Real>
  class Sequence
  {
  public:
    using Vector = Eigen::Vector<Real, Eigen::Dynamic>; /**< Vector of real numbers. */
    using Interval = struct Interval {
      Real    a; /**< Lower bound of the interval. */
      Real    b; /**< Upper bound of the interval. */
      Integer va; /**< Sign of the polynomial at the lower bound. */
      Integer vb; /**< Sign of the polynomial at the upper bound. */
      bool    a_on_root; /**< True if the lower bound is a root. */
      bool    b_on_root; /**< True if the upper bound is a root. */
    }; /**< Interval structure. */
    constexpr static const Real EPSILON{std::numeric_limits<Real>::epsilon()}; /**< Machine epsilon. */

  private:
    std::vector<Poly<Real>> m_sequence;  /**< Sturm sequence. */
    std::vector<Interval>   m_intervals; /**< Computed intervals. */
    Real                    m_a{0.0};    /**< Lower bound of the interval containing the roots. */
    Real                    m_b{0.0};    /**< Upper bound of the interval containing the roots. */

  public:
    /**
    * Class constructor for the Sturm sequence.
    */
    Sequence() {}

    /**
    * Class constructor for the Sturm sequence given the polynomial \f$ p(x) \f$.
    * \param[in] p Polynomial.
    */
    Sequence(Poly<Real> const & p) {this->build(p);}

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
    void build(Poly<Real> const & p) {
      this->m_intervals.clear();
      Poly<Real> dp, q, r;
      p.derivative(dp);
      this->m_sequence.clear();
      this->m_sequence.reserve(p.order());
      this->m_sequence.emplace_back(p);  this->m_sequence.back().adjust_degree();
      this->m_sequence.emplace_back(dp); this->m_sequence.back().adjust_degree();
      Integer n_sequence{1};
      while (true) {
        Sturm::divide<Real>(this->m_sequence[n_sequence - 1], this->m_sequence[n_sequence], q, r);
        if (r.order() <= 0) {break;}
        this->m_sequence.emplace_back(-r);
        ++n_sequence;
      }
      // Divide by GCD
      for (Integer i{0}; i < n_sequence; ++i) {
        Sturm::divide<Real>(this->m_sequence[i], this->m_sequence.back(), q, r);
        q.normalize();
        this->m_sequence[i] = q;
      }
      this->m_sequence.back().set_scalar(1);
    }


    /**
    * Get the length of the stored Sturm sequence.
    * \return Length of the stored Sturm sequence.
    */
    Integer length() const {return static_cast<Integer>(this->m_sequence.size());}

    /**
    * Get the \f$ i \f$-th polynomial of the stored Sturm sequence.
    * \return The \f$ i \f$-th polynomial of the stored Sturm sequence.
    */
    Poly<Real> const & get(Integer i) const {return this->m_sequence[i];}

    /**
    * Compute the sign variations of the stored Sturm sequence at \f$ x \f$.
    * \return The number of sign variations.
    */
    Integer sign_variations(Real x, bool & on_root ) const {
      Integer sign_var{0};
      Integer last_sign{0};
      Integer n_poly{Integer(this->m_sequence.size())};
      Real    v{this->m_sequence[0].evaluate(x)};
      on_root = false;
      if (v > 0) {last_sign = 1;}
      else if (v < 0) {last_sign = -1;}
      else {
        on_root   = true;
        last_sign = 0;
      }
      for (Integer i{1}; i < n_poly; ++i) {
        v = this->m_sequence[i].evaluate(x);
        if (v > 0) {
          if (last_sign == -1) {++sign_var;}
          last_sign = 1;
        } else if (v < 0) {
          if (last_sign == 1) {++sign_var;}
          last_sign = -1;
        }
      }
      return sign_var;
    }

    /**
    * Given an interval \f$ [a, b] \f$ compute the subintervals containing a single root.
    * \return The numbers of intervals (roots) found.
    */
    Integer separate_roots(Real a_in, Real b_in) {
      this->m_intervals.clear();
      this->m_intervals.reserve(this->m_sequence.size());

      Interval I_0, I_1;
      this->m_a = I_0.a = a_in;
      this->m_b = I_0.b = b_in;

      I_0.va = this->sign_variations(I_0.a, I_0.a_on_root);
      I_0.vb = this->sign_variations(I_0.b, I_0.b_on_root);

      Integer n_roots{std::abs(I_0.va - I_0.vb)};

      if (n_roots <= 1) {
        if (n_roots == 1 && !I_0.a_on_root && !I_0.b_on_root) {
          this->m_intervals.push_back(I_0);
        }
        if (I_0.a_on_root) {
          I_1.a  = I_1.b  = I_0.a;
          I_1.va = I_1.vb = I_0.va;
          I_1.a_on_root = I_1.b_on_root = true;
          this->m_intervals.push_back(I_1);
        }
        if (I_0.b_on_root) {
          I_1.a  = I_1.b  = I_0.b;
          I_1.va = I_1.vb = I_0.vb;
          I_1.a_on_root = I_1.b_on_root = true;
          this->m_intervals.push_back(I_1);
        }
        return static_cast<Integer>(m_intervals.size());
      }

      // Search intervals
      std::vector<Interval> I_stack; I_stack.clear(); I_stack.reserve(this->m_sequence.size());
      I_stack.push_back(I_0);
      while (I_stack.size() > 0) {
        I_0 = I_stack.back(); I_stack.pop_back();
        // Check if the interval has a single root
        n_roots = std::abs(I_0.va - I_0.vb);
        if (n_roots <= 1) {
          if (I_0.a_on_root) {
            I_0.b         = I_0.a;
            I_0.vb        = I_0.va;
            I_0.b_on_root = true;
            this->m_intervals.push_back(I_0);
          } else if (I_0.b_on_root) {
            I_0.a         = I_0.b;
            I_0.va        = I_0.vb;
            I_0.a_on_root = true;
            this->m_intervals.push_back(I_0);
          } else if (n_roots == 1) {
            this->m_intervals.push_back(I_0);
          }
        } else if (abs(I_0.b-I_0.a) <= static_cast<Real>(10.0)*EPSILON*std::max(static_cast<Real>(1.0),
          std::max(abs(I_0.b),abs(I_0.a)))) {
          I_1.a  = I_1.b  = I_0.a;
          I_1.va = I_1.vb = 0;
          I_1.a_on_root = I_1.b_on_root = true;
          I_stack.push_back(I_1);
        } else {
          Real    c{(I_0.a + I_0.b)/static_cast<Real>(2.0)};
          bool    c_on_root;
          Integer vc{this->sign_variations(c, c_on_root)};
          // Check the interval [a, c]
          if (I_0.va != vc || c_on_root || I_0.a_on_root) {
            if (c < I_1.b) { // Check if it is a true reduction
              I_1.a = I_0.a; I_1.va = I_0.va; I_1.a_on_root = I_0.a_on_root;
              I_1.b = c;    I_1.vb = vc;    I_1.b_on_root = c_on_root;
              I_stack.push_back(I_1);
            } else if (c_on_root) {
              I_1.a         = I_1.b         = c;
              I_1.a_on_root = I_1.b_on_root = true;
              I_1.va        = I_1.vb        = 0;
              I_stack.push_back(I_1);
            } else if (I_0.a_on_root) {
              I_1.a         = I_1.b         = I_0.a;
              I_1.a_on_root = I_1.b_on_root = true;
              I_1.va        = I_1.vb        = 0;
              I_stack.push_back(I_1);
            }
          }
          // Check the interval [c, b]
          if (I_0.vb != vc || I_0.b_on_root) {
            if (c > I_0.a) {
              I_1.a = c;    I_1.va = vc;    I_1.a_on_root = c_on_root;
              I_1.b = I_0.b; I_1.vb = I_0.vb; I_1.b_on_root = I_0.b_on_root;
              I_stack.push_back(I_1);
            } else if (I_0.b_on_root) {
              I_1.a         = I_1.b         = I_0.b;
              I_1.a_on_root = I_1.b_on_root = true;
              I_1.va        = I_1.vb        = 0;
              I_stack.push_back(I_1);
            }
          }
        }
      }
      // Sort intervals
      std::sort(this->m_intervals.begin(), this->m_intervals.end(),
        [](Interval const & I_a, Interval const & I_b) {return I_a.a < I_b.a;});
      return static_cast<Integer>(m_intervals.size());
    }

    /**
    * Compute an interval \f$ [a, b] \f$ that contains all the real roots and compute the subintervals
    * containing a single root.
    * \return The numbers of intervals (roots) found.
    */
    Integer separate_roots() {
      // Cauchy's bounds for roots
      Real leading_coeff{this->m_sequence[0].leading_coeff()};
      Real bound{static_cast<Real>(1.0) + this->m_sequence[0].cwiseAbs().maxCoeff() / std::abs(leading_coeff)};
      return separate_roots(-bound, bound);
    }

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

  }; // class Sequence

  /**
  * Print the Sturm sequence on an output stream.
  * \param[in] os Output stream.
  * \param[in] s Sturm sequence.
  * \return The output stream.
  * \tparam Real Scalar number type.
  */
  template <typename Real>
  inline std::ostream & operator<< (std::ostream & os, Sequence<Real> const & s) {
    // Print the Sturm sequence
    os << "Sturm sequence" << std::endl;
    for (Integer i{0}; i < s.length(); ++i) {
      os << "P_" << i << "(x) = " << s.get(i) << std::endl;
    }
    // Print the roots
    Integer n{s.roots_number()};
    if (n > 0) {
      os << "roots separation for interval [" << s.a() << "," << s.b() << "]" << std::endl;
      for (Integer i{0}; i < n; ++i) {
        typename Sequence<Real>::Interval const & I = s.interval(i);
        os << "I = [" << I.a << ", " << I.b << "], V = [" << I.va << ", " << I.vb << "]" << std::endl;
      }
    }
    return os;
  }

} // namespace Sturm

#endif // INCLUDE_STURM_SEQUENCE_HXX