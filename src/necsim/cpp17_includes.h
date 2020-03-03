// This file is part of necsim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @file cpp17_includes.h
 * @brief Contains headers for different systems, either using or not using boost and c++17 features.
 *
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */

#ifndef NECSIM_CPP17_INCLUDES_H
#define NECSIM_CPP17_INCLUDES_H




#ifdef WIN_INSTALL
#include <windows.h>
//#define sleep Sleep
#endif

// Provide support for C++ 14
#if __cplusplus < 201703L
#ifdef USING_BOOST

#include <boost/filesystem.hpp>

namespace fs = boost::filesystem;

#else

#include "ghc/filesystem.hpp"

namespace fs = ghc::filesystem;

#endif

#else // Support for C++ 17

#include <filesystem>

namespace fs = std::filesystem;

#endif // __cplusplus < 201703L
#endif //NECSIM_CPP17_INCLUDES_H
