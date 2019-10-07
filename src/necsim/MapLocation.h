//
// Created by sam on 07/09/19.
//

#ifndef NECSIM_MAPLOCATION_H
#define NECSIM_MAPLOCATION_H
#include <iostream>
namespace necsim
{
    struct MapLocation
    {
        long x;
        long y;
        long xwrap;
        long ywrap;

        MapLocation() : x(0), y(0), xwrap(0), ywrap(0)
        { }

        MapLocation(long x, long y, long xwrap, long ywrap) : x(x), y(y), xwrap(xwrap), ywrap(ywrap)
        { }

        /**
         * @brief Checks if the location is directly on the grid without wrapping (i.e. xwrap and ywrap are 0).
         * @return true if on the grid
         */
        bool isOnGrid() const;

        /**
         * @brief Equality operator for MapLocation
         * @param m the MapLocation object to compare against
         * @return true if the x, y, xwrap and ywrap are identical
         */
        bool operator==(MapLocation const &m) const;

        /**
         * @brief Inequality operator for MapLocation
         * @param m the MapLocation object to compare against
         * @return true if locations are not identical
         */
        bool operator!=(MapLocation const &m) const;

        /**
         * @brief Output operator for MapLocation
         * @param os the output stream to write to
         * @param m the object to write out
         * @return the output stream
         */
        friend std::ostream &operator<<(std::ostream &os, const MapLocation&m);

        /**
         * @brief Input operator for MapLocation
         * @param is the input stream to read from
         * @param m the object to read to
         * @return the input stream
         */
        friend std::istream &operator>>(std::istream &is, MapLocation &m);
    };
}
#endif //NECSIM_MAPLOCATION_H
