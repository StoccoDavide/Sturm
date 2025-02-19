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
  */
  class Poly : public Vector {

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
    * \param[in] a Scalar value.
    * \return Reference to the polynomial.
    */
    Poly & set_scalar(Real a);

    /**
    * Set the polynomial to a monomial with a given coefficient.
    * \param[in] a Coefficient of the monomial.
    * \return Reference to the polynomial.
    */
    Poly & set_monomial(Real a);

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
    std::string to_string() const;

    /**
    * Evaluate the polynomial at a given point.
    * \param[in] x Point at which to evaluate the polynomial.
    * \return Value of the polynomial at the given point.
    */
    Real evaluate(Real x) const;

    /**
    * Evaluate the derivative of the polynomial at a given point.
    * \param[in] x Point at which to evaluate the derivative.
    * \return Value of the derivative at the given point.
    */
    Real evaluate_derivative(Real x) const;

    /**
    * Evaluate the polynomial and its derivative at a given point.
    * \param[in] x Point at which to evaluate.
    * \param[out] p Value of the polynomial at the given point.
    * \param[out] Dp Value of the derivative at the given point.
    */
    void evaluate(Real x, Real & p, Real & Dp) const;

    /**
    * Get the leading coefficient of the polynomial.
    * \return Leading coefficient.
    */
    inline Real leading_coeff() const {return this->coeff(m_order-1);}

    /**
    * Compute the derivative of the polynomial.
    * \param[out] res Polynomial representing the derivative.
    */
    void derivative(Poly & res) const;

    /**
    * Compute the integral of the polynomial.
    * \param[out] res Polynomial representing the integral.
    */
    void integral(Poly & res) const;

    /**
    * Compute the integral of the polynomial with a given constant.
    * \param[out] res Polynomial representing the integral.
    * \param[in] c Constant value for the integral.
    */
    void integral(Poly & res, Real c) const;

    /**
    * Normalize the polynomial such that the maximum absolute coefficient is 1.
    * \return Scaling value used for normalization.
    */
    Real normalize();

    /**
    * Purge (set to 0) coefficients that are less than or equal to a given epsilon.
    * \param[in] epsi Epsilon value for purging coefficients.
    */
    void purge( Real epsi );

    /**
    * Adjust the degree of the polynomial to ensure the leading coefficient is non-zero.
    */
    void adjust_degree();

    /**
    * Count the number of sign variations in the polynomial coefficients.
    * \return Number of sign variations.
    */
    Integer sign_variations() const;

    /**
    * Change the polynomial in such a way \f$ p(x) = x^n + \sum_{i=0}^{n-1} a_i x^i \f$.
    */
    void make_monic();

    /**
    * Define the assignment with another polynomial.
    * \param[in] p Polynomial to assign from.
    * \return Reference to the polynomial.
    */
    Poly & operator= (Poly const & p);

    /**
    * Define the negation of the polynomial.
    * \return Negated polynomial.
    */
    Poly operator-();

    /**
    * Define the addition with another polynomial.
    * \param[in] p Polynomial to add.
    * \return Reference to the polynomial.
    */
    Poly & operator+= (Poly const & p);

    /**
    * Define the subtraction with another polynomial.
    * \param[in] p Polynomial to subtract.
    * \return Reference to the polynomial.
    */
    Poly & operator-= (Poly const & p);

    /**
    * Define the multiplication with another polynomial.
    * \param[in] p Polynomial to multiply.
    * \return Reference to the polynomial.
    */
    Poly & operator*= (Poly const & p);

    /**
    * Define the addition with a scalar.
    * \param[in] s Scalar to add.
    * \return Reference to the polynomial.
    */
    Poly & operator+= (Real s);

    /**
    * Define the subtraction with a scalar.
    * \param[in] s Scalar to subtract.
    * \return Reference to the polynomial.
    */
    Poly & operator-= (Real s);

    /**
    * Define the multiplication with a scalar.
    * \param[in] s Scalar to multiply.
    * \return Reference to the polynomial.
    */
    Poly & operator*= (Real s);

  }; // class Poly

  /**
  * Define the sum between two polynomials.
  * \param[in] a First polynomial to sum.
  * \param[in] b Second polynomial to sum.
  * \return The operation result.
  */
  Poly operator+ (Poly const & a, Poly const & b);

  /**
  * Define the sum between a polynomial and a scalar.
  * \param[in] a Polynomial to sum.
  * \param[in] b Scalar to sum.
  * \return The operation result.
  */
  Poly operator+ (Poly const & a, Real b);

  /**
  * Define the sum between a polynomial and a scalar.
  * \param[in] a Scalar to sum.
  * \param[in] b Polynomial to sum.
  * \return The operation result.
  */
  Poly operator+ (Real a, Poly const & b);

  /**
  * Define the difference between two polynomials.
  * \param[in] a First polynomial to subtract.
  * \param[in] b Second polynomial to subtract.
  * \return The operation result.
  */
  Poly operator- (Poly const & a, Poly const & b);

  /**
  * Define the difference between a scalar and a polynomial.
  * \param[in] a Polynomial to subtract.
  * \param[in] b Scalar to subtract.
  * \return The operation result.
  */
  Poly operator- (Poly const & a, Real b);

  /**
  * Define the difference between a scalar and a polynomial.
  * \param[in] a Scalar to subtract.
  * \param[in] b Polynomial to subtract.
  * \return The operation result.
  */
  Poly operator- (Real a, Poly const & b);

  /**
  * Define the multiplication between two polynomials.
  * \param[in] a First polynomial to multiply.
  * \param[in] b Second polynomial to multiply.
  * \return The operation result.
  */
  Poly operator* (Poly const & a, Poly const & b);


  /**
  * Define the multiplication between a scalar and a polynomial.
  * \param[in] a Polynomial to multiply.
  * \param[in] b Scalar to multiply.
  * \return The operation result.
  */
  Poly operator* (Poly const & a, Real b);

  /**
  * Define the multiplication between a scalar and a polynomial.
  * \param[in] a Scalar to multiply.
  * \param[in] b Polynomial to multiply.
  * \return The operation result.
  */
  Poly operator* (Real a, Poly const & b);

  /**
  * Divide the polynomial \f$ a(x) \f$ by \f$ b(x) \f$ with remainder \f$ r(x) \f$ and quotient \f$ q(x) \f$.
  * \param[in] a Polynomial to divide.
  * \param[in] b Polynomial to divide by.
  * \param[out] q Quotient polynomial.
  * \param[out] r Remainder polynomial.
  */
  void divide(
    Poly const & a,
    Poly const & b,
    Poly       & q,
    Poly       & r
  );

  /**
  * Compute the greatest common divisor of two polynomials.
  * \param[in] a First polynomial.
  * \param[in] b Second polynomial.
  * \param[out] g Greatest common divisor polynomial.
  * \param[in] epsi Epsilon value for purging coefficients.
  */
  void GCD(Poly const & a, Poly const & b, Poly & g, Real epsi = EPSILON);

  /**
  * Print the polynomial on an output stream.
  * \param[in] os Output stream.
  * \param[in] p Polynomial.
  * \return The output stream.
  */
  std::ostream & operator<< (std::ostream & os, Poly const & p);

} // namespace Sturm

#endif // INCLUDE_STURM_POLY_HXX
