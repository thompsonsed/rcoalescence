//This file is part of necsim project which is released under MIT license.
//See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @author Samuel Thompson
 * @date 16/08/2017
 * @file ActivityMap.h
 * @brief Contains the ActivityMap, which inherits from Matrix and adds a few extra parameters.
 *
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */

#ifndef ACTIVITYMAP_H
#define ACTIVITYMAP_H

#include <cstring>
#include <string>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <memory>
#include "Map.h"
#include "RNGController.h"

/**
 * @brief Contains the routines for importing activity maps and getting a cell value from the map.
 *
 * Activity maps can be reproduction probability maps or death maps.
 */
class ActivityMap
{
protected:
    // Matrix containing the relative activity probabilities
    Map<double> activity_map;
    // Path to the map file
    string map_file;
    // Maximum value across the map
    double max_val;
    // If true,
    bool null_map;
    // The fine map offsets and the sample map dimensions
    unsigned long offset_x, offset_y, x_dim, y_dim;
    // Random number generator
    shared_ptr<RNGController> random;

    // Function pointer for our reproduction map checker
    typedef bool (ActivityMap::*rep_ptr)(const unsigned long &x, const unsigned long &y, const long &xwrap,
                                         const long &ywrap);

    // once setup will contain the end check function to use for this simulation.
    rep_ptr activity_map_checker_fptr;
public:
    ActivityMap() : activity_map(), offset_x(0), offset_y(0), x_dim(0), y_dim(0), random(make_shared<RNGController>()),
                    activity_map_checker_fptr(nullptr)
    {
        map_file = "none";
        max_val = 0;
        null_map = true;
    }

    /**
     * @brief Gets if the map is a null reproduction map, with no matrix storage in memory necessary.
     * @return true if the map is "null" or "none"
     */
    bool isNull();

    /**
     * @brief Imports the map file from the given path
     * Requires the dimensions to be identical as the fine map file dimensions
     * @param file_name the path to the reproduction map to import
     * @param size_x the x dimensions of the map file
     * @param size_y the y dimensions of the map file
     */
    void import(string file_name, unsigned long size_x, unsigned long size_y, shared_ptr<RNGController> random_in);

    /**
     * @brief Correctly sets the reproduction function to either rejectionSampleNull or rejectionSample
     * depending on if a reproduction map is used or not.
     */
    void setActivityFunction();

    /**
     * @brief Sets the offsets for the reproduction map from the sample grid.
     * @param x_offset the x offset from the sample grid
     * @param y_offset the y offset from the sample grid
     * @param xdim the x dimension of the sample grid
     * @param ydim the y dimension of the sample grid
     */
    void setOffsets(const unsigned long &x_offset, const unsigned long &y_offset,
                    const unsigned long &xdim, const unsigned long &ydim);

    /**
     * @brief Returns true for all cell values
     * Function to be pointed to in cases where there is no reproduction map
     * @param random_number random number object to pass forward
     * @param x x coordinate of the lineage on the sample grid
     * @param y y coordinate of the lineage on the sample grid
     * @param xwrap x wrapping of the lineage
     * @param ywrap y wrapping of the lineage
     * @return true always
     */
    bool rejectionSampleNull(const unsigned long &x, const unsigned long &y, const long &xwrap, const long &ywrap);

    /**
     * @brief Returns true for all cell values
     * Function to be pointed to in cases where there is no reproduction map
     * @param random_number random number object to pass forward
     * @param x x coordinate of the lineage on the sample grid
     * @param y y coordinate of the lineage on the sample grid
     * @param xwrap x wrapping of the lineage
     * @param ywrap y wrapping of the lineage
     * @return true always
     */
    bool rejectionSample(const unsigned long &x, const unsigned long &y, const long &xwrap, const long &ywrap);

    /**
     * @brief Gets the value of the reproduction map at that location
     * @param x x coordinate of the lineage on the sample grid
     * @param y y coordinate of the lineage on the sample grid
     * @param xwrap x wrapping of the lineage
     * @param ywrap y wrapping of the lineage
     * @return value of the reproduction map at the required location
     */
    double getVal(const unsigned long &x, const unsigned long &y, const long &xwrap, const long &ywrap);

    /**
     * @brief Tests if the random action occurs
     * @param x x coordinate of the lineage on the sample grid
     * @param y y coordinate of the lineage on the sample grid
     * @param xwrap x wrapping of the lineage
     * @param ywrap y wrapping of the lineage
     * @return value of the reproduction map at the required location
     */
    bool actionOccurs(const unsigned long &x, const unsigned long &y, const long &xwrap, const long &ywrap);

    /**
     * @brief Standardises probability values from 0-1.
     */
    void standardiseValues();

    /**
     * @brief Get the value at the specified index.
     * @param rows the row index to obtain
     * @param cols the column index to obtain
     * @return the value at the specified row and column
     */
    double get(const unsigned long &rows, const unsigned long &cols)
    {
        return activity_map.get(rows, cols);
    }

    /**
     * @brief Equality operator
     * @param rm the ActivityMap object to copy from
     * @return the self ActivityMap object
     */
    ActivityMap &operator=(const ActivityMap &rm)
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

    /**
     * @brief Operator for outputting to an ostream.
     * @param os the ostream to output to
     * @param r the ActivityMap to read from
     * @return the os object
     */
    friend ostream &operator<<(ostream &os, ActivityMap &r)
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

    /**
     * @brief Operator for inputting from an istream.
     * @param is the istream to input from
     * @param r the ActivityMap to input to
     * @return the is object
     */
    friend istream &operator>>(istream &is, ActivityMap &r)
    {
        is.ignore();
        getline(is, r.map_file);
        unsigned long col, row;
        is >> col >> row;
        is >> r.offset_x >> r.offset_y >> r.x_dim >> r.y_dim;
        r.import(r.map_file, col, row, shared_ptr<RNGController>());
        return is;
    }

};

#endif //ACTIVITYMAP_H
