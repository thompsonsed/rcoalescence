//
// Created by sam on 07/09/19.
//

#include "MapLocation.h"
namespace necsim
{
    bool MapLocation::isOnGrid() const
    {
        return xwrap == 0 && ywrap == 0;
    }

    bool MapLocation::operator==(MapLocation const &m) const
    {
        return x == m.x && y == m.y && xwrap == m.xwrap && ywrap == m.ywrap;
    }

    bool MapLocation::operator!=(MapLocation const &m) const
    {
        return !(this->operator==(m));
    }

    std::ostream &operator<<(std::ostream &os, const MapLocation &m)
    {
        os << m.x << "," << m.y << "," << m.xwrap << "," << m.ywrap << std::endl;
        return os;
    }

    std::istream &operator>>(std::istream &is, MapLocation &m)
    {
        char delim;
        is >> m.x >> delim >> m.y >> delim >> m.xwrap >> delim >> m.ywrap;
        return is;
    }
}
