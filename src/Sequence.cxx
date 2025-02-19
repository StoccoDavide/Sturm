/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sturm project is distributed under the BSD 2-Clause License.                              *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * davide.stocco@unitn.it                                             enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "Sturm.hh"

namespace Sturm {

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void Sequence::build(Poly const & p) {
    this->m_intervals.clear();
    Poly dp, q, r;
    p.derivative(dp);
    this->m_sequence.clear();
    this->m_sequence.reserve(p.order());
    this->m_sequence.emplace_back(p);  this->m_sequence.back().adjust_degree();
    this->m_sequence.emplace_back(dp); this->m_sequence.back().adjust_degree();
    Integer n_sequence{1};
    while (true) {
      Sturm::divide(this->m_sequence[n_sequence - 1], this->m_sequence[n_sequence], q, r);
      if (r.order() <= 0) {break;}
      this->m_sequence.emplace_back(-r);
      ++n_sequence;
    }
    // Divide by GCD
    for (Integer i{0}; i < n_sequence; ++i) {
      Sturm::divide(this->m_sequence[i], this->m_sequence.back(), q, r);
      q.normalize();
      this->m_sequence[i] = q;
    }
    this->m_sequence.back().set_scalar(1);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Integer Sequence::sign_variations(Real x, bool & on_root ) const {
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

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


  Integer Sequence::separate_roots(Real a_in, Real b_in) {
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
      } else if (abs(I_0.b-I_0.a) <= 10*EPSILON*std::max(1.0, std::max(abs(I_0.b),abs(I_0.a)))) {
        I_1.a  = I_1.b  = I_0.a;
        I_1.va = I_1.vb = 0;
        I_1.a_on_root = I_1.b_on_root = true;
        I_stack.push_back(I_1);
      } else {
        Real    c{(I_0.a + I_0.b)/2.0};
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

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Integer Sequence::separate_roots() {
    // Cauchy's bounds for roots
    Real leading_coeff{this->m_sequence[0].leading_coeff()};
    Real bound{1.0 + this->m_sequence[0].cwiseAbs().maxCoeff() / std::abs(leading_coeff)};
    return separate_roots(-bound, bound);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void
  Sequence::refine_roots() {
    this->m_function = std::function<void(Real x, Real & out)>([this](Real x, Real & out) {
      out = this->m_sequence[0].evaluate(x);
    });
    this->m_roots.resize(this->m_intervals.size());
    Integer n{0};
    for (auto & I : this->m_intervals) {
      Real & r{this->m_roots.coeffRef(n++)};
      if (I.a_on_root) {r = I.a;}
      else if (I.b_on_root) {r = I.b;}
      else {
        this->m_solver.bounds(I.a, I.b);
        this->m_solver.solve(this->m_function, 0.0, r);
        STURM_ASSERT_WARNING(this->m_solver.converged(),
          "Sturm::Sequence::refine_roots(...): failed at interval n = " << n);
      }
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  std::ostream & operator<< (std::ostream & os, Sequence const & s) {
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
        Sequence::Interval const & I = s.interval(i);
        os << "I = [" << I.a << ", " << I.b << "], V = [" << I.va << ", " << I.vb << "]" << std::endl;
      }
    }
    return os;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

} // namespace Sturm
