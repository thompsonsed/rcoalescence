// This file is part of rcoalescence project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.
/**
 * @author Sam Thompson
 * @file RLogging.cpp
 * @brief Controls logging from R.
 *
 * Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */
#include <Rcpp/r/headers.h>
#include <Rcpp.h>

#include "RLogging.h"

using namespace std;

void writeInfo(string message)
{
    if(logging_mode)
    {
        Rcpp::Rcout << message << flush;
    }
}

void writeWarning(string message)
{
    writeInfo(message);
}

void writeError(string message)
{
    writeInfo(message);
}

void writeCritical(string message)
{
    writeInfo(message);
}


