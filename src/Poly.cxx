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

namespace Sturm
{

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  std::string Poly::to_string() const {

    if (this->order() <= 0) {return "(empty polynomial)";}
    if (this->order() == 1) {return std::to_string(this->coeff(0));};
    if ((*this).cwiseAbs().maxCoeff() == 0) {return "0";};

    bool empty{true}; // true means that all the coefficients are zero
    std::string s{""}; // sign
    Real        c{0}; // coefficient
    std::string e{""}; // exponent
    std::string res{""};

    // Check if the first coefficient is different from zero
    if (this->coeff(0) != 0) {
      res = std::to_string(this->coeff(0));
      empty = false;
    }

    for (Integer i{1}; i < this->order(); ++i) {
      // If the coefficient is negative...
      if (this->coeff(i) < 0) {
        // and if the previous coefficients are null...
        if (empty) {
          s = ""; // ...do not write the sign
          c = this->coeff(i);
          empty = false;
        } else {
          s = " - "; // ...otherwise write the sign as a separator
          c = -this->coeff(i); // ...and change the sign of the coefficient
        }

      } else if (this->coeff(i) > 0) {
        // If the coefficient is positive...
        c = this->coeff(i);
        // and if the previous coefficients are null...
        if (empty) {
          s = ""; // ...do not write the sign
          empty = false;
        } else {
          s = " + "; // ...otherwise write the sign as a separator
        }

      } else {
        // If the coefficient is zero, skip to the next
        continue;
      }

      // If the exponent is 1, do not write it
      if (i == 1) {e = "x";}
      else {e = "x^" + std::to_string(i);}

      // If the coefficient is 1, do not write it
      if (c == 1) {res += s; res += e;}
      else {res += s + std::to_string(c) + e;}
    }
    return res;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Poly & Poly::set_scalar(Real a) {
    this->resize(1);
    this->coeffRef(0) = a;
    this->m_order     = 1;
    return *this;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Poly & Poly::set_monomial(Real a) {
    this->resize(2);
    this->coeffRef(0) = a;
    this->coeffRef(1) = 1;
    this->m_order     = 2;
    return *this;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Real Poly::evaluate(Real x) const {
    Integer n{this->m_order - 1};
    Real res{this->coeff(n)};
    while (n-- > 0) {res = res*x+this->coeff(n);}
    return res;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Real Poly::evaluate_derivative(Real x) const {
    Integer n{this->m_order - 1};
    Real Dp{this->coeff(n)*n};
    while (--n > 0) {Dp = Dp*x+this->coeff(n)*n;}
    return Dp;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void Poly::evaluate(Real x, Real & p, Real & Dp) const {
    Integer n{this->m_order - 1};
    p  = this->coeff(n);
    Dp = this->coeff(n)*n;
    while (--n > 0) {
      p  = p*x  + this->coeff(n);
      Dp = Dp*x + this->coeff(n)*n;
    }
    p = p*x+this->coeff(0);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void Poly::derivative(Poly & res) const {
    res.resize(this->m_order - 1);
    for (Integer i{1}; i < this->m_order; ++i) {
      res.coeffRef(i - 1) = i * this->coeff(i);
    }
    res.m_order = this->m_order - 1;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void Poly::integral(Poly & res) const {
    res.resize(this->m_order + 1);
    res.coeffRef(0) = 0;
    for (Integer i{1}; i <= this->m_order; ++i) {
      res.coeffRef(i) = this->coeff(i - 1)/i;
    }
    res.m_order = this->m_order + 1;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Real Poly::normalize() {
    Real scale{this->cwiseAbs().maxCoeff()};
    if (scale > 0.0) {this->to_eigen() /= scale;}
    return scale;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void Poly::purge(Real eps) {
    if (this->m_order > 0) {
      Real max{this->cwiseAbs().maxCoeff()};
      if (max < 1) {max = 1;}
      Real eps_tmp{eps*max};
      for (Integer i{0}; i < this->m_order; ++i) {
        Real & c_i{this->coeffRef(i)};
        if (std::abs(c_i) <= eps_tmp) {c_i = 0;}
      }
    }
    this->adjust_degree();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void Poly::adjust_degree() {
    while (this->m_order > 0 && this->coeff(this->m_order - 1) == 0) {--this->m_order;}
    this->conservativeResize(this->m_order);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Integer Poly::sign_variations() const {
    Integer sign_var{0};
    Integer last_sign{0};
    for (Integer i{0}; i < this->m_order; ++i) {
      Real v = this->coeff(i);
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

  void Poly::make_monic() {
    this->to_eigen() /= this->coeff(this->m_order - 1);
    this->coeffRef(this->m_order - 1) = 1;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Poly & Poly::operator= (Poly const & b) {
    this->resize(b.m_order);
    this->to_eigen().noalias() = b.to_eigen();
    this->m_order = b.m_order;
    return *this;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Poly Poly::operator- () {
    return Poly(-this->to_eigen());
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Poly & Poly::operator+= (Poly const & b) {
    Integer max_order{std::max(this->m_order, b.m_order)};
    Integer min_order{std::min(this->m_order, b.m_order)};
    this->conservativeResize(max_order); // Resize the vector without destroying the content
    this->head( min_order ).noalias() += b.head(min_order);
    Integer n_tail{b.m_order - this->m_order};
    if (n_tail > 0) {this->tail(n_tail).noalias() = b.tail(n_tail);}
    this->m_order = max_order;
    return *this;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Poly & Poly::operator+= (Real b) {
    if (this->m_order > 0) this->coeffRef(0) += b;
    else {
      this->resize(1);
      this->coeffRef(0) = b;
      this->m_order = 1;
    }
    return *this;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Poly & Poly::operator-= (Poly const & b) {
    Integer max_order{std::max(this->m_order, b.m_order)};
    Integer min_order{std::min(this->m_order, b.m_order)};
    this->conservativeResize(max_order); // Resize the vector without destroying the content
    this->head( min_order ).noalias() -= b.head(min_order);
    Integer n_tail{b.m_order - this->m_order};
    if (n_tail > 0) {this->tail(n_tail).noalias() = -b.tail(n_tail);}
    this->m_order = max_order;
    return *this;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Poly & Poly::operator-= (Real b) {
    if (this->m_order > 0) {
      this->coeffRef(0) -= b;
    } else {
      this->resize(1);
      this->coeffRef(0) = -b;
      this->m_order     = 1;
    }
    return *this;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Poly & Poly::operator*= (Poly const & b) {
    Vector a(this->to_eigen());
    Integer new_order{this->m_order + b.m_order - 1};
    this->resize(new_order);
    this->setZero();
    for(Integer i{0}; i<this->m_order; ++i) {
      for(Integer j{0}; j<b.m_order; ++j) {
        this->coeffRef(i+j) += a.coeff(i) * b.coeff(j);
    }}
    this->m_order = new_order;
    return *this;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Poly & Poly::operator*= (Real b) {
    this->to_eigen() *= b;
    return *this;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Poly operator+ (Poly const & a, Poly const & b) {
    Integer max_order{std::max(a.order(), b.order())};
    Integer min_order{std::min(a.order(), b.order())};
    Poly res(max_order);
    res.head( min_order ).noalias() = a.head(min_order) + b.head(min_order);
    Integer n_tail{max_order - min_order};
    if (n_tail > 0) {
      if (a.order() > b.order()) {res.tail(n_tail).noalias() = a.tail(n_tail);}
      else {res.tail(n_tail).noalias() = b.tail(n_tail);}
    }
    return res;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Poly operator+ (Poly const & a, Real b) {
    Integer max_order{std::max(a.order(), 1)};
    Poly res(max_order);
    if (a.order() > 0) {
      res.coeffRef(0) = a.coeff(0) + b;
      if (a.order() > 1) res.tail( a.order() - 1).noalias() = a.tail(a.order() - 1);
    } else {
      res.coeffRef(0) = b;
    }
    return res;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Poly operator+ (Real a, Poly const & b) {
    Integer max_order{std::max(b.order(), 1)};
    Poly res(max_order);
    if (b.order() > 0) {
      res.coeffRef(0) = a + b.coeff(0);
      if (b.order() > 1) res.tail(b.order() - 1).noalias() = b.tail(b.order() - 1);
    } else {
      res.coeffRef(0) = a;
    }
    return res;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Poly operator- (Poly const & a, Poly const & b) {
    Integer max_order{std::max(a.order(), b.order())};
    Integer min_order{std::min(a.order(), b.order())};
    Poly res(max_order);
    res.head( min_order ).noalias() = a.head(min_order) - b.head(min_order);
    Integer n_tail{max_order - min_order};
    if (n_tail > 0) {
      if (a.order() > b.order()) {res.tail(n_tail).noalias() = a.tail(n_tail);}
      else {res.tail(n_tail).noalias() = -b.tail(n_tail);}
    }
    return res;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Poly operator- (Poly const & a, Real b) {
    Integer max_order{std::max(a.order(), 1)};
    Poly res(max_order);
    if (a.order() > 0) {
      res.coeffRef(0) = a.coeff(0) - b;
      if (a.order() > 1) res.tail( a.order() - 1).noalias() = a.tail(a.order() - 1);
    } else {
      res.coeffRef(0) = -b;
    }
    return res;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Poly operator- (Real a, Poly const & b) {
    Integer max_order{std::max(b.order(), 1)};
    Poly res(max_order);
    if (b.order() > 0) {
      res.coeffRef(0) = a - b.coeff(0);
      if (b.order() > 1) {res.tail(b.order() - 1).noalias() = -b.tail(b.order() - 1);}
    } else {
      res.coeffRef(0) = a;
    }
    return res;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Poly operator* (Poly const & a, Poly const & b) {
    Poly res(a.order() + b.order() - 1);
    for(Integer i{0}; i < a.order(); ++i) {
      for(Integer j{0}; j < b.order(); ++j) {
        res.coeffRef(i+j) += a.coeff(i) * b.coeff(j);
    }}
    return res;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Poly operator* (Real a, Poly const & b) {
    Poly res(b.order());
    res.noalias() = a*b.to_eigen();
    return res;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Poly operator* (Poly const & a, Real b) {
    Poly res(a.order());
    res.noalias() = a.to_eigen()*b;
    return res;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void divide(Poly const & a, Poly const & b, Poly & q, Poly & r) {
    // Scale polynomials as a_norm(x) = a(x) / scale_a, and b_norm(x) = b(x) / scale_b
    Poly a_norm(a), b_norm(b);
    Real scale_a = a_norm.normalize();
    Real scale_b = b_norm.normalize();

    // a_norm(x) = b_norm(x) * q(x) + r(x)
    r = a_norm;
    Real leading_b_norm{b_norm.leading_coeff()};
    Integer d{r.order() - b_norm.order()};
    if (d < 0) {
      // a_norm(x) = b_norm(x) + r(x)
      q.set_scalar(1);
      r = a_norm - b_norm;
    } else {
      Integer r_degree{r.degree()};
      q.set_order(d+1);

      STURM_ASSERT(std::abs(leading_b_norm) > EPSILON,
        "Sturm::Poly::divide(...): leading coefficient of b(x) is 0.");

      while (d >= 0 && r_degree >= 0) {
        Real leading_r = {r(r_degree)};
        Real ratio{leading_r/leading_b_norm};
        q.coeffRef(d) = ratio;
        r.segment(d, b_norm.degree()).noalias() -= ratio*b_norm.head(b_norm.degree());
        r.coeffRef(r_degree) = 0;
        --r_degree;
        --d;
      }

      // Do not purge remainder: this can be done externally with r.purge(eps);
      r.adjust_degree();
    }

    /*
    Scale back polinomials:
    - a_norm(x) = b_norm(x) * q(x) + r(x)
    - a(x) / scale_a = b(x) / scale_b * q(x) + r(x)
    - a(x) = b(x) * (scale_a/scale_b) * q(x) + scale_a*r(x)
    */
    q *= scale_a/scale_b;
    r *= scale_a;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void GCD(Poly const & a, Poly const & b, Poly & g, Real eps)
  {
    if (b.order() > 0) {
      Poly q, r;
      Sturm::divide(a, b, q, r);
      r.purge(eps);
      Sturm::GCD(b, r, g, eps);
    } else {
      g = a;
    }
    g.normalize();
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  std::ostream & operator<< (std::ostream & os, Poly const & p) {
    os << p.to_string();
    return os;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

}