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

#include "ActivityMap.h"
#include <utility>

bool ActivityMap::isNull()
{
    return null_map;
}

void ActivityMap::import(string file_name, unsigned long size_x, unsigned long size_y, shared_ptr<RNGController> random_in)
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

void ActivityMap::setOffsets(const unsigned long &x_offset, const unsigned long &y_offset,
                             const unsigned long &xdim, const unsigned long &ydim)
{
    offset_x = x_offset;
    offset_y = y_offset;
    x_dim = xdim;
    y_dim = ydim;
}

bool
ActivityMap::rejectionSampleNull(const unsigned long &x, const unsigned long &y, const long &xwrap, const long &ywrap)
{
    return true;
}

bool ActivityMap::rejectionSample(const unsigned long &x, const unsigned long &y, const long &xwrap, const long &ywrap)
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






