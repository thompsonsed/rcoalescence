// This file is part of necsim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @author Samuel Thompson
 * @file DoubleComparison.h
 *
 * @copyright <a href="https://opensource.org/licenses/MIT">MIT Licence.</a>
 * @brief Simple comparison between floating points.
 *
 * Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
 */

#ifndef NECSIM_DOUBLECOMPARISON_H
#define NECSIM_DOUBLECOMPARISON_H

/**
 * @brief Compares two doubles and returns a boolean of whether they are equal, within the epsilon deviation.
 * This is useful for floating point errors in saving and reading doubles from file.
 * @param d1 the first double.
 * @param d2 the second double.
 * @param epsilon the deviation within which the values are assumed to be equal.
 * @return true if the doubles are within epsilon of each other
 */
bool doubleCompare(double d1, double d2, double epsilon);

/**
 * @brief Compares two doubles and returns a boolean of whether they are equal, within the epsilon deviation.
 * This is useful for floating point errors in saving and reading doubles from file.
 * Overloaded version for long doubles.
 * @param d1 the first double.
 * @param d2 the second double.
 * @param epsilon the deviation within which the values are assumed to be equal.
 * @return true if the doubles are within epsilon of each other
 */
bool doubleCompare(long double d1, long double d2, long double epsilon);

/**
 * @brief Compares two doubles and returns a boolean of whether they are equal, within the epsilon deviation.
 * This is useful for floating point errors in saving and reading doubles from file.
 * Overloaded version for long doubles and double epsilon.
 * @param d1 the first double.
 * @param d2 the second double.
 * @param epsilon the deviation within which the values are assumed to be equal.
 * @return true if the doubles are within epsilon of each other
 */
bool doubleCompare(long double d1, long double d2, double epsilon);

#endif //NECSIM_DOUBLECOMPARISON_H
