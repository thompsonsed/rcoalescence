// This file is part of necsim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.


/**
 * @author Samuel Thompson
 * @file DataMask.h
 * @brief  Contains DataMask for describing the spatial sampling pattern on a landscape.
 *
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */
#ifndef SPECIATIONCOUNTER_DataMask_H
#define SPECIATIONCOUNTER_DataMask_H

// Forward declaration of Landscape class
class Landscape;

#include <string>
#include <memory>

#include "SimParameters.h"
#include "Map.h"


// Class which contains the DataMask object, telling us where to sample from within the habitat map.
/**
 * @brief Contains the DataMask object, a Matrix of booleans describing the spatial sampling pattern.
 */
class DataMask
{
protected:
    // the file to read in from
    string inputfile;
    // True if the sample mask is true everywhere
    bool isNullSample;
    // True if the grid is smaller than the sample mask
    bool isGridOffset;
    unsigned long x_offset, y_offset;
    // Stores the size of the grid which is stored as a full species species_id_list
    unsigned long x_dim, y_dim;
    // Stores the size of the samplemask from which spatially sampling is read
    unsigned long mask_x_dim, mask_y_dim;

    // Function pointer for obtaining the proportional sampling from the sample mask.
    typedef double (DataMask::*fptr)(const long &x, const long &y, const long &xwrap, const long &ywrap) const;

    fptr getProportionfptr;
public:
    /** A binary grid telling whether or not the cell should be sampled.*/
    Map<bool> sample_mask;
    /** Exact grid for determining sampling proportion.*/
    Map<double> sample_mask_exact;

    /**
     * @brief The DataMask constructor.
     */
    DataMask();

    ~DataMask() = default;

    /**
     * @brief Returns if the simulation is using the a null samplemask, and therefore does not need to store the full
     * sample grid in memory.
     * @return true if using a null samplemask
     */
    bool isNull();

    /**
     * @brief Sets the parameters for the datamask, including the dimensions of the map, the offsets from the grid and
     * the dimensions of the grid itself for recalculating coordinates.
     * @param sim_parameters Simulation parameter to set the data mask from
     */
    void setup(shared_ptr<SimParameters> sim_parameters);

    bool checkCanUseDefault(shared_ptr<SimParameters> sim_parameters);

    /**
     * @brief Imports the sample mask as a boolean mask and sets the relevant sample mask dimensions.
     * Should only be called if the import is actually required (i.e. the map is not null or none).
     * @param xdim the x dimension of the grid area
     * @param ydim the y dimension of the grid area
     * @param mask_xdim the x dimension of the sample map file
     * @param mask_ydim the y dimension of the sample map file
     * @param xoffset the x offset of the grid area from the sample map file
     * @param yoffset the y offset of the grid area from the sample map file
     * @param inputfile_in the path to the sample map file
     */
    void importBooleanMask(unsigned long xdim, unsigned long ydim, unsigned long mask_xdim, unsigned long mask_ydim,
                           unsigned long xoffset, unsigned long yoffset, string inputfile_in);

    /**
     * @brief Imports the boolean map object.
     */
    void doImport();

    /**
     * @brief Sets the map dimensions and the getVal function pointer.
     */
    void completeBoolImport();

    /**
     * @brief Sets up the null sampling map.
     */
    void setupNull(shared_ptr<SimParameters> mapvarin);

    /**
     * @brief Imports the specified file for the sampling percentage within each cell.
     *
     * The map should consist of floating points representing the relative sampling rates in each cell. Note that the
     * actual sampling proportion is equal to the cell value multiplied by global deme sampling proportion.
     * @param mapvarin the SimParameters object containing the samplemask file location and dimensions
     */
    void importSampleMask(shared_ptr<SimParameters> mapvarin);

    /**
     * @brief Calculates the matrix value at the provided x, y location.
     * If everywhere is sampled, simply returns true, as no sample_mask will be stored in memory.
     * This is to save RAM where possible.
     *
     * @param x the x position on the grid
     * @param y the y position on the grid
     * @param xwrap the number of x wraps
     * @param ywrap the number of y wraps
     * @return the sample_mask value at x,y (or true if the file was "null")
     */
    bool getVal(const long &x, const long &y, const long &xwrap, const long &ywrap) const;

    /**
     * @brief Separate return function for always returning 1.0 as density value.
     * @param x the x position on the grid
     * @param y the y position on the grid
     * @param xwrap the number of x wraps around the map
     * @param ywrap the number of y wraps around the map
     * @return the sample_mask_exact value at (x, y)
     */
    double getNullProportion(const long &x, const long &y, const long &xwrap, const long &ywrap) const;

    /**
     * @brief Returns the exact value from the spatial sampling map, for calculating the proportion of individuals
     * to be sampled in each cell.
     * @note this function assumes that the file is not "null" and the exact sampling mask has been imported. No error
     * checks on these conditions are performed except in debugging mode.
     * @param x the x position on the grid
     * @param y the y position on the grid
     * @param xwrap the number of x wraps around the map
     * @param ywrap the number of y wraps around the map
     * @return the sample_mask_exact value at (x, y)
     */
    double getBoolProportion(const long &x, const long &y, const long &xwrap, const long &ywrap) const;

    /**
     * @brief Returns the exact value from the spatial sampling map, for calculating the proportion of individuals
     * to be sampled in each cell.
     * @note this function assumes that the file is not "null" and the exact sampling mask has been imported. No error
     * checks on these conditions are performed except in debugging mode.
     * @param x the x position on the grid
     * @param y the y position on the grid
     * @param xwrap the number of x wraps around the map
     * @param ywrap the number of y wraps around the map
     * @return the sample_mask_exact value at (x, y)
     */
    double getSampleProportion(const long &x, const long &y, const long &xwrap, const long &ywrap) const;

    /**
     * @brief Returns the exact value from the spatial sampling map, as returned by the pointer function.
     * @param x the x position on the grid
     * @param y the y position on the grid
     * @param xwrap the number of x wraps around the map
     * @param ywrap the number of y wraps around the map
     * @return the sample_mask_exact value at (x, y)
     */
    double getExactValue(const long &x, const long &y, const long &xwrap, const long &ywrap) const;

    /**
     * @brief Converts the spatial map into the boolean grid required for continued simulation.
     * This is done so that the faster boolean accesses are possible.
     * @param map1 the map object to obtain density values from
     * @param deme_sampling the proportion of individuals to sample
     * @param generation the generation individuals are added at
     */
    void convertBoolean(shared_ptr<Landscape> map1, const double &deme_sampling, const double &generation);

    /**
     * @brief Removes the spatial mask from memory. This should be performed if no more map expansions are required.
     */
    void clearSpatialMask();

    /**
     * @brief Converts the coordinates back into the grid format. Changes the values in the provided variables to
     * be correct.
     *
     * @param x the x value to convert
     * @param y the y value to convert
     * @param x_wrap the xwrap variable to place the value into
     * @param y_wrap the ywrap variable to place the value into
     */
    void recalculateCoordinates(long &x, long &y, long &x_wrap, long &y_wrap);
};

#endif //SPECIATIONCOUNTER_DataMask_H
