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

#ifndef INCLUDE_STURM_HH
#define INCLUDE_STURM_HH

// C++ standard libraries
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

// Eigen library
#include <Eigen/Dense>

// Print Sturm errors
#ifndef STURM_ERROR
#define STURM_ERROR(MSG)                \
  {                                     \
    std::ostringstream os;              \
    os << MSG;                          \
    throw std::runtime_error(os.str()); \
  }
#endif

// Assert for Sturm
#ifndef STURM_ASSERT
#define STURM_ASSERT(COND, MSG) \
  if (!(COND))                  \
  {                             \
    STURM_ERROR(MSG);           \
  }
#endif

// Warning for Sturm
#ifndef STURM_WARNING
#define STURM_WARNING(MSG)         \
  {                                \
    std::cout << MSG << std::endl; \
  }
#endif

// Warning assert for Sturm
#ifndef STURM_ASSERT_WARNING
#define STURM_ASSERT_WARNING(COND, MSG) \
  if (!(COND))                          \
  {                                     \
    STURM_WARNING(MSG);                 \
  }
#endif

/**
* \brief Namespace for the Sturm library.
*/
namespace Sturm
{

  using Real    = double; /**< Real number type. */
  using Integer = int;    /**< Integer number type. */
  using Vector  = Eigen::Vector<Real, Eigen::Dynamic>; /**< Vector of real numbers. */

  static const Real EPSILON = std::numeric_limits<Real>::epsilon(); /**< Machine epsilon. */

  /**
  * Print Sturm library information on a string.
  * \return A string with the Sturm library information.
  */
  std::string Info() {
    std::ostringstream os;
    os
      << "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *" << std::endl
      << "* Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *" << std::endl
      << "*                                                                                               *" << std::endl
      << "* The Sturm project is distributed under the BSD 2-Clause License.                              *" << std::endl
      << "*                                                                                               *" << std::endl
      << "* Davide Stocco                                                               Enrico Bertolazzi *" << std::endl
      << "* University of Trento                                                     University of Trento *" << std::endl
      << "* davide.stocco@unitn.it                                             enrico.bertolazzi@unitn.it *" << std::endl
      << "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *" << std::endl;
    return os.str();
  }

  /**
  * Print Sturm library information on a stream.
  * \param[in] os Output stream.
  */
  void Info(std::ostream & os) {os << Info();}

} // namespace Sturm

#include "Sturm/Poly.hxx"
#include "Sturm/Sequence.hxx"

#endif // INCLUDE_STURM_HH
