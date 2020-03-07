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
// #ifdef R_COMPILER
// #undef Realloc
// #ifndef R_Realloc
// #define R_Realloc(p,n,t) (t *) R_chk_realloc( (void *)(p), (size_t)((n) * sizeof(t)) )
// #endif
// #include <shlobj.h>
// #undef Free
// #endif
#include <windows.h>
#include <winerror.h>
#ifndef ERROR_FILE_TOO_LARGE
#define ERROR_FILE_TOO_LARGE 223L
#endif

//#define sleep Sleep
#endif

// Provide support for C++ 14
#if __cplusplus < 201703L
#ifdef USING_BOOST

#include <boost/filesystem.hpp>

namespace fs = boost::filesystem;

#else
#ifndef ERROR_FILE_TOO_LARGE
#define ERROR_FILE_NOT_FOUND 2L
#define ERROR_PATH_NOT_FOUND 3L
#define ERROR_ACCESS_DENIED 5L
#define ERROR_NO_MORE_FILES 18L
#define ERROR_NOT_SUPPORTED 50L
#define ERROR_INVALID_PARAMETER 87L
#define ERROR_CALL_NOT_IMPLEMENTED 120L
#define ERROR_INVALID_NAME 123L
#define ERROR_ALREADY_EXISTS 183L
#define ERROR_FILE_TOO_LARGE 223L
#define ERROR_PRIVILEGE_NOT_HELD 1314L
#endif

#include <ghc/filesystem.hpp>

namespace fs = ghc::filesystem;

#endif

#else // Support for C++ 17

#include <filesystem>

namespace fs = std::filesystem;

#endif // __cplusplus < 201703L
#endif //NECSIM_CPP17_INCLUDES_H