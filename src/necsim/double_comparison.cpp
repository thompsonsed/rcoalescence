// This file is part of necsim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @author Samuel Thompson
 * @file DoubleComparison.cpp
 *
 * @copyright <a href="https://opensource.org/licenses/MIT">MIT Licence.</a>
 * @brief Simple comparison between floating points.
 *
 * Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
 */
#include <cmath>

using namespace std;

bool doubleCompare(double d1, double d2, double epsilon)
{
    return (abs(float(d1 - d2)) < epsilon);
}

bool doubleCompare(long double d1, long double d2, long double epsilon)
{
    return abs((d1 - d2)) < epsilon;
}

bool doubleCompare(long double d1, long double d2, double epsilon)
{
    return abs((d1 - d2)) < epsilon;
}
