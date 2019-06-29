// This file is part of necsim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @author Samuel Thompson
 * @file Cell.h
 * @brief Basic container for location data.
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */

#ifndef CELL_H
#define CELL_H

/**
 * @class Cell
 * @brief Simple structure containing the x and y positions of a cell
 */
struct Cell
{
    long x;
    long y;

    /**
     * @brief Overloading equality operator
     * @param c the Cell containing the values to overload
     * @return the cell with the new values
     */
    Cell &operator=(Cell const &c)
    = default;

    /**
     * @brief Equality operator for Cell
     * @param c the Cell object to compare against
     * @return true if the x and y locations are identical
     */
    bool operator==(Cell const &c)
    {
        return x == c.x && y == c.y;
    }

    /**
     * @brief Inequality operator for Cell
     * @param c the Cell object to compare against
     * @return true if the x and y locations are not identical
     */
    bool operator!=(Cell const &c)
    {
        return !(this->operator==(c));
    }
};

/**
 * @brief Calculates the distance between two cells
 *
 * @param c1 Cell containing one point
 * @param c2 Cell containing second point
 * @return the distance between the two points
 */
double distanceBetweenCells(Cell &c1, Cell &c2);

#endif // CELL_H
