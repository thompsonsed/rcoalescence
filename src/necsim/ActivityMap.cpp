//This file is part of necsim project which is released under MIT license.
//See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @author Samuel Thompson
 * @date 16/08/2017
 * @file ActivityMap.cpp
 * @brief Contains the ActivityMap, which inherits from Matrix and adds a few extra parameters.
 *
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */

#include <utility>
#include <numeric>

#include "ActivityMap.h"

namespace necsim
{
    bool ActivityMap::isNull()
    {
        return null_map;
    }

    void ActivityMap::import(string file_name,
                             unsigned long size_x,
                             unsigned long size_y,
                             shared_ptr<RNGController> random_in)
    {
        random = std::move(random_in);
        map_file = file_name;
        if(file_name == "null" || file_name == "none")
        {
            null_map = true;
        }
        else
        {
            null_map = false;
            activity_map.setSize(size_y, size_x);
            activity_map.import(file_name);
            for(unsigned long y = 0; y < activity_map.getRows(); y++)
            {
                for(unsigned long x = 0; x < activity_map.getCols(); x++)
                {
                    if(activity_map.get(y, x) > max_val)
                    {
                        max_val = activity_map.get(y, x);
                    }
                }
            }
            activity_map.close();
        }
        setActivityFunction();
    }

    void ActivityMap::setActivityFunction()
    {
        if(null_map)
        {
            activity_map_checker_fptr = &ActivityMap::rejectionSampleNull;
        }
        else
        {
            activity_map_checker_fptr = &ActivityMap::rejectionSample;
        }
    }

    void ActivityMap::setOffsets(const unsigned long &x_offset,
                                 const unsigned long &y_offset,
                                 const unsigned long &xdim,
                                 const unsigned long &ydim)
    {
        offset_x = x_offset;
        offset_y = y_offset;
        x_dim = xdim;
        y_dim = ydim;
    }

    bool ActivityMap::rejectionSampleNull(const unsigned long &x,
                                          const unsigned long &y,
                                          const long &xwrap,
                                          const long &ywrap)
    {
        return true;
    }

    bool ActivityMap::rejectionSample(const unsigned long &x,
                                      const unsigned long &y,
                                      const long &xwrap,
                                      const long &ywrap)
    {
        return random->d01() <= getVal(x, y, xwrap, ywrap);
    }

    double ActivityMap::getVal(const unsigned long &x, const unsigned long &y, const long &xwrap, const long &ywrap)
    {
        unsigned long x_ref = x + (xwrap * x_dim) + offset_x;
        unsigned long y_ref = y + (ywrap * y_dim) + offset_y;
        return activity_map.get(y_ref, x_ref);
    }

    bool ActivityMap::actionOccurs(const unsigned long &x, const unsigned long &y, const long &xwrap, const long &ywrap)
    {
        return (this->*activity_map_checker_fptr)(x, y, xwrap, ywrap);
    }

    void ActivityMap::standardiseValues()
    {
        if(!isNull())
        {
            double max_value = 0;
            for(unsigned long i = 0; i < activity_map.getRows(); i++)
            {
                for(unsigned long j = 0; j < activity_map.getCols(); j++)
                {
                    if(activity_map.get(i, j) > max_value)
                    {
                        max_value = activity_map.get(i, j);
                    }
                }
            }
            if(max_value == 0)
            {
                throw FatalException("Activity map does not contain any probability values.");
            }
            activity_map /= max_value;
        }
    }

    double ActivityMap::get(const unsigned long &rows, const unsigned long &cols)
    {
        return activity_map.get(rows, cols);
    }

    double ActivityMap::getMean() const
    {
        return std::accumulate(activity_map.begin(), activity_map.end(), 0.0);
    }

    ActivityMap &ActivityMap::operator=(const ActivityMap &rm)
    {
        this->activity_map = rm.activity_map;
        this->map_file = rm.map_file;
        this->max_val = rm.max_val;
        this->null_map = rm.null_map;
        this->offset_x = rm.offset_x;
        this->offset_y = rm.offset_y;
        this->x_dim = rm.x_dim;
        this->y_dim = rm.y_dim;
        this->activity_map_checker_fptr = rm.activity_map_checker_fptr;
        return *this;
    }

    ostream &operator<<(ostream &os, ActivityMap &r)
    {
        os << r.map_file << "\n";
        os << r.activity_map.getCols() << "\n";
        os << r.activity_map.getRows() << "\n";
        os << r.offset_x << "\n";
        os << r.offset_y << "\n";
        os << r.x_dim << "\n";
        os << r.y_dim << "\n";
        return os;
    }

    istream &operator>>(istream &is, ActivityMap &r)
    {
        is.ignore();
        getline(is, r.map_file);
        unsigned long col, row;
        is >> col >> row;
        is >> r.offset_x >> r.offset_y >> r.x_dim >> r.y_dim;
        r.import(r.map_file, col, row, shared_ptr<RNGController>());
        return is;
    }


}
