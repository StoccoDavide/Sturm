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

#ifndef INCLUDE_STURM_POLY_HXX
#define INCLUDE_STURM_POLY_HXX

namespace Sturm
{

  /**
  * \brief Polynomial class.
  *
  * This class implements a polynomial of the form \f$ p(x) = \sum_{i=0}^n a_i x^i \f$.
  * \tparam Real Scalar number type.
  */
  template <typename Real>
  class Poly : public Eigen::Vector<Real, Eigen::Dynamic>
  {
  public:
    using Vector = Eigen::Vector<Real, Eigen::Dynamic>; /**< Vector of real numbers. */
    constexpr static const Real EPSILON{std::numeric_limits<Real>::epsilon()}; /**< Machine epsilon. */

  private:
    Integer m_order; /**< Polynomial order. */

    /**
    * Access to the Eigen class.
    * \return The Eigen class.
    */
    Vector & to_eigen() {return *static_cast<Vector*>(this);}

  public:
    /**
    * Class constructor from a vector.
    */
    Poly() : m_order(0) {}

    /**
    * Class constructor from a vector.
    * \param[in] order Order of the polynomial.
    */
    explicit Poly(Integer order) : m_order(order) {
      this->resize(order);
      this->setZero();
    }

    /**
    * Class constructor from a vector.
    * \param[in] c Vector of coefficients.
    */
    Poly(Poly const & c) : Vector(c), m_order(c.m_order) {}

    /**
    * Class constructor from a vector.
    * \param[in] c Vector of coefficients.
    */
    explicit Poly(Vector const & c) {
      this->resize( c.size() );
      this->Vector::operator= (c); // to avoid recursion
      this->m_order = Integer(c.size());
    }

    /**
    * Access to the Eigen const class.
    * \return The Eigen const class.
    */
    Vector const & to_eigen() const {return *static_cast<Vector const *>(this);}

    /**
    * Set the order of the polynomial.
    * \param[in] order Order of the polynomial.
    */
    void
    set_order(Integer order) {
      this->m_order = order;
      this->resize( order );
      this->setZero();
    }

    /**
    * Set the degree of the polynomial.
    * \param[in] degree Degree of the polynomial.
    */
    void set_degree(Integer degree)
    {
      this->m_order = degree + 1;
      this->resize(this->m_order);
      this->setZero();
    }

    /**
    * Set the polynomial to a scalar value.
    * \param[in] s Scalar value.
    * \return Reference to the polynomial.
    */
    Poly & set_scalar(Real s) {
      this->resize(1);
      this->coeffRef(0) = s;
      this->m_order = 1;
      return *this;
    }

    /**
    * Set the polynomial to a monomial with a given coefficient.
    * \param[in] a Coefficient of the monomial.
    * \return Reference to the polynomial.
    */
    Poly & set_monomial(Real a) {
      this->resize(2);
      this->coeffRef(0) = a;
      this->coeffRef(1) = 1;
      this->m_order = 2;
      return *this;
    }

    /**
    * Get the coefficients of the polynomial.
    * \return Coefficients at the given index.
    */
    Vector coeffs() const {return this->to_eigen();}

    /**
    * Get the degree of the polynomial.
    * \return Degree of the polynomial.
    */
    Integer degree() const {return this->m_order - 1;}

    /**
    * Get the order of the polynomial.
    * \return Order of the polynomial.
    */
    Integer order() const {return this->m_order;}

    /**
    * Convert the polynomial to a string representation.
    * \return String representation of the polynomial.
    */
    std::string to_string() const {

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

    /**
    * Evaluate the polynomial at a given point.
    * \param[in] x Point at which to evaluate the polynomial.
    * \return Value of the polynomial at the given point.
    */
    Real evaluate(Real x) const {
      Integer n{this->m_order - 1};
      Real res{this->coeff(n)};
      while (n-- > 0) {res = res*x+this->coeff(n);}
      return res;
    }

    /**
    * Evaluate the derivative of the polynomial at a given point.
    * \param[in] x Point at which to evaluate the derivative.
    * \return Value of the derivative at the given point.
    */
    Real evaluate_derivative(Real x) const {
      Integer n{this->m_order - 1};
      Real Dp{this->coeff(n)*n};
      while (--n > 0) {Dp = Dp*x+this->coeff(n)*n;}
      return Dp;
    }

    /**
    * Evaluate the polynomial and its derivative at a given point.
    * \param[in] x Point at which to evaluate.
    * \param[out] p Value of the polynomial at the given point.
    * \param[out] Dp Value of the derivative at the given point.
    */
    void evaluate(Real x, Real & p, Real & Dp) const {
      Integer n{this->m_order - 1};
      p  = this->coeff(n);
      Dp = this->coeff(n)*n;
      while (--n > 0) {
        p  = p*x  + this->coeff(n);
        Dp = Dp*x + this->coeff(n)*n;
      }
      p = p*x+this->coeff(0);
    }

    /**
    * Get the leading coefficient of the polynomial.
    * \return Leading coefficient.
    */
    Real leading_coeff() const {return this->coeff(m_order-1);}

    /**
    * Compute the derivative of the polynomial.
    * \param[out] res Polynomial representing the derivative.
    */
    void derivative(Poly & res) const {
      res.resize(this->m_order - 1);
      for (Integer i{1}; i < this->m_order; ++i) {
        res.coeffRef(i - 1) = i * this->coeff(i);
      }
      res.m_order = this->m_order - 1;
    }

    /**
    * Compute the integral of the polynomial.
    * \param[out] res Polynomial representing the integral.
    */
    void integral(Poly & res) const {
      res.resize(this->m_order + 1);
      res.coeffRef(0) = 0;
      for (Integer i{1}; i <= this->m_order; ++i) {
        res.coeffRef(i) = this->coeff(i - 1)/i;
      }
      res.m_order = this->m_order + 1;
    }

    /**
    * Compute the integral of the polynomial with a given constant.
    * \param[out] res Polynomial representing the integral.
    * \param[in] c Constant value for the integral.
    */
    void integral(Poly & res, Real c) const {
      res.resize(this->m_order + 1);
      res.coeffRef(0) = c;
      for (Integer i{1}; i <= this->m_order; ++i) {
        res.coeffRef(i) = this->coeff(i - 1)/i;
      }
      res.m_order = this->m_order + 1;
    }

    /**
    * Normalize the polynomial such that the maximum absolute coefficient is 1.
    * \return Scaling value used for normalization.
    */
    Real normalize() {
      Real scale{this->cwiseAbs().maxCoeff()};
      if (scale > static_cast<Real>(0.0)) {this->to_eigen() /= scale;}
      return scale;
    }

    /**
    * Purge (set to 0) coefficients that are less than or equal to a given epsilon.
    * \param[in] eps Epsilon value for purging coefficients.
    */
    void purge(Real eps) {
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

    /**
    * Adjust the degree of the polynomial to ensure the leading coefficient is non-zero.
    */
    void adjust_degree() {
      while (this->m_order > 0 && this->coeff(this->m_order - 1) == 0) {--this->m_order;}
      this->conservativeResize(this->m_order);
    }

    /**
    * Count the number of sign variations in the polynomial coefficients.
    * \return Number of sign variations.
    */
    Integer sign_variations() const {
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

    /**
    * Change the polynomial in such a way \f$ p(x) = x^n + \sum_{i=0}^{n-1} a_i x^i \f$.
    */
    void make_monic() {
      this->to_eigen() /= this->coeff(this->m_order - 1);
      this->coeffRef(this->m_order - 1) = 1;
    }

    /**
    * Define the assignment with another polynomial.
    * \param[in] p Polynomial to assign from.
    * \return Reference to the polynomial.
    */
    Poly & operator= (Poly const & p) {
      this->resize(p.m_order);
      this->to_eigen().noalias() = p.to_eigen();
      this->m_order = p.m_order;
      return *this;
    }

    /**
    * Define the negation of the polynomial.
    * \return Negated polynomial.
    */
    Poly operator- () {
      return Poly(-this->to_eigen());
    }

    /**
    * Define the addition with another polynomial.
    * \param[in] p Polynomial to add.
    * \return Reference to the polynomial.
    */
    Poly & operator+= (Poly const & p) {
      Integer max_order{std::max(this->m_order, p.m_order)};
      Integer min_order{std::min(this->m_order, p.m_order)};
      this->conservativeResize(max_order); // Resize the vector without destroying the content
      this->head( min_order ).noalias() += p.head(min_order);
      Integer n_tail{p.m_order - this->m_order};
      if (n_tail > 0) {this->tail(n_tail).noalias() = p.tail(n_tail);}
      this->m_order = max_order;
      return *this;
    }

    /**
    * Define the subtraction with another polynomial.
    * \param[in] p Polynomial to subtract.
    * \return Reference to the polynomial.
    */
    Poly & operator-= (Poly const & p) {
      Integer max_order{std::max(this->m_order, p.m_order)};
      Integer min_order{std::min(this->m_order, p.m_order)};
      this->conservativeResize(max_order); // Resize the vector without destroying the content
      this->head( min_order ).noalias() -= p.head(min_order);
      Integer n_tail{p.m_order - this->m_order};
      if (n_tail > 0) {this->tail(n_tail).noalias() = -p.tail(n_tail);}
      this->m_order = max_order;
      return *this;
    }

    /**
    * Define the multiplication with another polynomial.
    * \param[in] p Polynomial to multiply.
    * \return Reference to the polynomial.
    */
    Poly & operator*= (Poly const & p) {
      Vector a(this->to_eigen());
      Integer new_order{this->m_order + p.m_order - 1};
      this->resize(new_order);
      this->setZero();
      for(Integer i{0}; i<this->m_order; ++i) {
        for(Integer j{0}; j < p.m_order; ++j) {
          this->coeffRef(i+j) += a.coeff(i) * p.coeff(j);
      }}
      this->m_order = new_order;
      return *this;
    }

    /**
    * Define the addition with a scalar.
    * \param[in] s Scalar to add.
    * \return Reference to the polynomial.
    */
    Poly & operator+= (Real s) {
      if (this->m_order > 0) this->coeffRef(0) += s;
      else {
        this->resize(1);
        this->coeffRef(0) = s;
        this->m_order = 1;
      }
      return *this;
    }

    /**
    * Define the subtraction with a scalar.
    * \param[in] s Scalar to subtract.
    * \return Reference to the polynomial.
    */
    Poly & operator-= (Real s) {
      if (this->m_order > 0) {
        this->coeffRef(0) -= s;
      } else {
        this->resize(1);
        this->coeffRef(0) = -s;
        this->m_order = 1;
      }
      return *this;
    }

    /**
    * Define the multiplication with a scalar.
    * \param[in] s Scalar to multiply.
    * \return Reference to the polynomial.
    */
    Poly & operator*= (Real s) {
      this->to_eigen() *= s;
      return *this;
    }

  }; // class Poly

  /**
  * Define the sum between two polynomials.
  * \param[in] p_1 First polynomial to sum.
  * \param[in] p_2 Second polynomial to sum.
  * \return The operation result.
  * \tparam Real Scalar number type.
  */
  template <typename Real>
  inline Poly<Real> operator+ (Poly<Real> const & p_1, Poly<Real> const & p_2) {
    Integer max_order{std::max(p_1.order(), p_2.order())};
    Integer min_order{std::min(p_1.order(), p_2.order())};
    Poly<Real> res(max_order);
    res.head( min_order ).noalias() = p_1.head(min_order) + p_2.head(min_order);
    Integer n_tail{max_order - min_order};
    if (n_tail > 0) {
      if (p_1.order() > p_2.order()) {res.tail(n_tail).noalias() = p_1.tail(n_tail);}
      else {res.tail(n_tail).noalias() = p_2.tail(n_tail);}
    }
    return res;
  }

  /**
  * Define the sum between a polynomial and a scalar.
  * \param[in] p Polynomial to sum.
  * \param[in] s Scalar to sum.
  * \return The operation result.
  * \tparam Real Scalar number type.
  */
  template <typename Real>
  inline Poly<Real> operator+ (Poly<Real> const & p, Real s) {
    Integer max_order{std::max(p.order(), 1)};
    Poly<Real> res(max_order);
    if (p.order() > 0) {
      res.coeffRef(0) = p.coeff(0) + s;
      if (p.order() > 1) res.tail(p.order() - 1).noalias() = p.tail(p.order() - 1);
    } else {
      res.coeffRef(0) = s;
    }
    return res;
  }

  /**
  * Define the sum between a polynomial and a scalar.
  * \param[in] s Scalar to sum.
  * \param[in] p Polynomial to sum.
  * \return The operation result.
  * \tparam Real Scalar number type.
  */
  template <typename Real>
  inline Poly<Real> operator+ (Real s, Poly<Real> const & p) {
    Integer max_order{std::max(p.order(), 1)};
    Poly<Real> res(max_order);
    if (p.order() > 0) {
      res.coeffRef(0) = s + p.coeff(0);
      if (p.order() > 1) res.tail(p.order() - 1).noalias() = p.tail(p.order() - 1);
    } else {
      res.coeffRef(0) = s;
    }
    return res;
  }

  /**
  * Define the difference between two polynomials.
  * \param[in] p_1 First polynomial to subtract.
  * \param[in] p_2 Second polynomial to subtract.
  * \return The operation result.
  * \tparam Real Scalar number type.
  */
  template <typename Real>
  inline Poly<Real> operator- (Poly<Real> const & p_1, Poly<Real> const & p_2) {
    Integer max_order{std::max(p_1.order(), p_2.order())};
    Integer min_order{std::min(p_1.order(), p_2.order())};
    Poly<Real> res(max_order);
    res.head( min_order ).noalias() = p_1.head(min_order) - p_2.head(min_order);
    Integer n_tail{max_order - min_order};
    if (n_tail > 0) {
      if (p_1.order() > p_2.order()) {res.tail(n_tail).noalias() = p_1.tail(n_tail);}
      else {res.tail(n_tail).noalias() = -p_2.tail(n_tail);}
    }
    return res;
  }

  /**
  * Define the difference between a scalar and a polynomial.
  * \param[in] p Polynomial to subtract.
  * \param[in] s Scalar to subtract.
  * \return The operation result.
  * \tparam Real Scalar number type.
  */
  template <typename Real>
  inline Poly<Real> operator- (Poly<Real> const & p, Real s) {
    Integer max_order{std::max(p.order(), 1)};
    Poly<Real> res(max_order);
    if (p.order() > 0) {
      res.coeffRef(0) = p.coeff(0) - s;
      if (p.order() > 1) res.tail( p.order() - 1).noalias() = p.tail(p.order() - 1);
    } else {
      res.coeffRef(0) = -s;
    }
    return res;
  }

  /**
  * Define the difference between a scalar and a polynomial.
  * \param[in] s Scalar to subtract.
  * \param[in] p Polynomial to subtract.
  * \return The operation result.
  * \tparam Real Scalar number type.
  */
  template <typename Real>
  inline Poly<Real> operator- (Real s, Poly<Real> const & p) {
    Integer max_order{std::max(p.order(), 1)};
    Poly<Real> res(max_order);
    if (p.order() > 0) {
      res.coeffRef(0) = s - p.coeff(0);
      if (p.order() > 1) {res.tail(p.order() - 1).noalias() = -p.tail(p.order() - 1);}
    } else {
      res.coeffRef(0) = s;
    }
    return res;
  }

  /**
  * Define the multiplication between two polynomials.
  * \param[in] p_1 First polynomial to multiply.
  * \param[in] p_2 Second polynomial to multiply.
  * \return The operation result.
  * \tparam Real Scalar number type.
  */
  template <typename Real>
  inline Poly<Real> operator* (Poly<Real> const & p_1, Poly<Real> const & p_2) {
    Poly<Real> res(p_1.order() + p_2.order() - 1);
    for(Integer i{0}; i < p_1.order(); ++i) {
      for(Integer j{0}; j < p_2.order(); ++j) {
        res.coeffRef(i+j) += p_1.coeff(i) * p_2.coeff(j);
    }}
    return res;
  }

  /**
  * Define the multiplication between a scalar and a polynomial.
  * \param[in] p Polynomial to multiply.
  * \param[in] s Scalar to multiply.
  * \return The operation result.
  * \tparam Real Scalar number type.
  */
  template <typename Real>
  inline Poly<Real> operator* (Real p, Poly<Real> const & s) {
    Poly<Real> res(s.order());
    res.noalias() = p*s.to_eigen();
    return res;
  }

  /**
  * Define the multiplication between a scalar and a polynomial.
  * \param[in] s Scalar to multiply.
  * \param[in] p Polynomial to multiply.
  * \return The operation result.
  * \tparam Real Scalar number type.
  */
  template <typename Real>
  inline Poly<Real> operator* (Poly<Real> const & s, Real p) {
    Poly<Real> res(s.order());
    res.noalias() = s.to_eigen()*p;
    return res;
  }

  /**
  * Divide the polynomial \f$ p_1(x) \f$ by \f$ p_2(x) \f$ with remainder \f$ r(x) \f$ and quotient \f$ q(x) \f$.
  * \param[in] p_1 Polynomial to divide.
  * \param[in] p_2 Polynomial to divide by.
  * \param[out] q Quotient polynomial.
  * \param[out] r Remainder polynomial.
  * \tparam Real Scalar number type.
  */
  template <typename Real>
  inline void divide(Poly<Real> const & p_1, Poly<Real> const & p_2, Poly<Real> & q, Poly<Real> & r) {
    // Scale polynomials as p_1_norm(x) = p_1(x) / scale_p_1, and p_2_norm(x) = p_2(x) / scale_p_2
    Poly<Real> p_1_norm(p_1), p_2_norm(p_2);
    Real scale_p_1{p_1_norm.normalize()};
    Real scale_p_2{p_2_norm.normalize()};

    // p_1_norm(x) = p_2_norm(x) * q(x) + r(x)
    r = p_1_norm;
    Real leading_b_norm{p_2_norm.leading_coeff()};
    Integer d{r.order() - p_2_norm.order()};
    if (d < 0) {
      // p_1_norm(x) = p_2_norm(x) + r(x)
      q.set_scalar(1);
      r = p_1_norm - p_2_norm;
    } else {
      Integer r_degree{r.degree()};
      q.set_order(d+1);

      STURM_ASSERT(std::abs(leading_b_norm) > Poly<Real>::EPSILON,
        "Sturm::Poly::divide(...): leading coefficient of p_2(x) is 0.");

      while (d >= 0 && r_degree >= 0) {
        Real leading_r = {r(r_degree)};
        Real ratio{leading_r/leading_b_norm};
        q.coeffRef(d) = ratio;
        r.segment(d, p_2_norm.degree()).noalias() -= ratio*p_2_norm.head(p_2_norm.degree());
        r.coeffRef(r_degree) = 0;
        --r_degree;
        --d;
      }

      // Do not purge remainder: this can be done externally with r.purge(eps);
      r.adjust_degree();
    }

    /*
    Scale back polinomials:
    - p_1_norm(x) = p_2_norm(x) * q(x) + r(x)
    - p_1(x) / scale_p_1 = p_2(x) / scale_p_2 * q(x) + r(x)
    - p_1(x) = p_2(x) * (scale_p_1/scale_p_2) * q(x) + scale_p_1*r(x)
    */
    q *= scale_p_1/scale_p_2;
    r *= scale_p_1;
  }

  /**
  * Compute the greatest common divisor of two polynomials.
  * \param[in] p_1 First polynomial.
  * \param[in] p_2 Second polynomial.
  * \param[out] gcd Greatest common divisor polynomial.
  * \param[in] eps Epsilon value for purging coefficients.
  * \tparam Real Scalar number type.
  */
  template <typename Real>
  inline void GCD(Poly<Real> const & p_1, Poly<Real> const & p_2, Poly<Real> & gcd, Real eps = Poly<Real>::EPSILON)
  {
    if (p_2.order() > 0) {
      Poly<Real> q, r;
      Sturm::divide(p_1, p_2, q, r);
      r.purge(eps);
      Sturm::GCD(p_2, r, gcd, eps);
    } else {
      gcd = p_1;
    }
    gcd.normalize();
  }

  /**
  * Print the polynomial on an output stream.
  * \param[in] os Output stream.
  * \param[in] p Polynomial.
  * \return The output stream.
  * \tparam Real Scalar number type.
  */
  template <typename Real>
  inline std::ostream & operator<< (std::ostream & os, Poly<Real> const & p) {
    os << p.to_string();
    return os;
  }

} // namespace Sturm

#endif // INCLUDE_STURM_POLY_HXX
