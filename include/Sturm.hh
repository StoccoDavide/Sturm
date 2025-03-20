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

// Default integer type
#ifndef STURM_DEFAULT_INTEGER_TYPE
#define STURM_DEFAULT_INTEGER_TYPE int
#endif

/**
* \brief Namespace for the Sturm library.
*
* This namespace contains all the classes and functions of the Sturm library for the computation of the
* Sturm sequences and the greatest common divisor polynomials. Other basic functions to manipulate
* polynomials and sequences are also included.
*/
namespace Sturm {

  /**
  * \brief The Integer type as used for the API.
  *
  * The Integer type, \c \#define the preprocessor symbol \c STURM_DEFAULT_INTEGER_TYPE. The default
  * value is \c int.
  */
  using Integer = STURM_DEFAULT_INTEGER_TYPE;

}

#endif // INCLUDE_STURM_HH
