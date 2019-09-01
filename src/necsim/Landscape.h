//This file is part of necsim project which is released under MIT license.
//See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details
/**
 * @author Samuel Thompson
 * @date 31/08/16
 * @file Landscape.h
 *
 * @brief Contains the Landscape object for easy referencing of the respective coarse and fine map within the same coordinate system.
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */

#ifndef LANDSCAPE_H
#define LANDSCAPE_H

#include <string>
#include <cstdio>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include <memory>
#include <boost/filesystem.hpp>

#include "Map.h"
#include "DataMask.h"
#include "SimParameters.h"

using namespace std;

/**
 * @brief Imports the provided given file to the provided matrix, given the x and y dimensions.
 * The integer value in the final matrix is multiplied by the scalar to generate the final matrix. Note that this
 * doubles the memory usage.
 * @param map_file the path to the map file to import
 * @param map_in the matrix to change the values of
 * @param map_x the x dimension of the matrix
 * @param map_y the y dimension of the matrix
 * @param scalar the scalar to multiply all values in the final matrix by (before rounding to integer)
 * @return the maximum value from the imported matrix
 */
uint32_t importToMapAndRound(string map_file, Map<uint32_t> &map_in, unsigned long map_x,
                             unsigned long map_y, double scalar);

/**
 * @brief Gets the x coordinate of the archimedes spiral
 * @param centre_x the x coordinate of the spiral centre
 * @param centre_y the y coordinate of the spiral centre
 * @param radius the radius of the spiral
 * @param theta the theta of rotation of the spiral
 */
unsigned long archimedesSpiralX(const double &centre_x, const double &centre_y, const double &radius,
                                const double &theta);

/**
 * @brief Gets the y coordinate of the archimedes spiral
 * @param centre_x the x coordinate of the spiral centre
 * @param centre_y the y coordinate of the spiral centre
 * @param radius the radius of the spiral
 * @param theta the theta of rotation of the spiral
 */
unsigned long archimedesSpiralY(const double &centre_x, const double &centre_y, const double &radius,
                                const double &theta);

/**
 * @brief Calculates the distance between two points.
 * @param start_x the start x coordinate
 * @param start_y the start y coordinate
 * @param end_x the end x coordinate
 * @param end_y the end y coodinate
 * @return the distance between the points
 */
double calculateDistance(const double &start_x, const double &start_y, const double &end_x, const double &end_y);

/**
 * @class Landscape
 * @brief Contains all maps and provides the functions for accessing a grid cell in the correct temporal and spacial
 * location.
 *
 * @details The function runDispersal() also provides the move routine, provided two alternative methods for moving
 * individuals. Contains routines for easy setting up and switching between the different coordinate systems required.
 * Set the map parameters with setDims(), import the map files with calcFineMap(), calcCoarseMap() etc, then set up
 * the landscape type using setLandscape() and sethistorical().
 * Usage is then by runDispersal() for running a dispersal kernel on the landscape, and then getVal() to obtain the
 * density at the desired coordinates. All coordinates should be given in reference to the simulation grid, and offsets
 * for the fine and coarse map are calculated automatically.
 *
 */
class Landscape
{
protected:
    // The map files which are read in (or generated if running with "null" as the map file".
    // historical maps are meant for before any deforestation occured, whereas the other maps are intended for modern day maps.
    // A linear transformation from modern to historical maps is used, approaching the habitat_change_rate variable times the difference between the historical and modern maps.
    // Once the gen_since_historical number of generations has been reached, the map will jump to the historical condition.
    // the finer grid for the area around the sample area.
    Map<uint32_t> fine_map;
    // the historical finer map.
    Map<uint32_t> historical_fine_map;
    // the coarser grid for the wider zone.
    Map<uint32_t> coarse_map;
    // the historical coarser map.
    Map<uint32_t> historical_coarse_map;
    // for importing and storing the simulation set-up options.
    shared_ptr<SimParameters> mapvars;
    // the minimum values for each dimension for offsetting.
    long fine_x_min{}, fine_y_min{}, coarse_x_min{}, coarse_y_min{};
    // the maximum values for each dimension for offsetting.
    long fine_x_max{}, fine_y_max{}, coarse_x_max{}, coarse_y_max{};
    // the offsetting of the map in FINE map units.
    long fine_x_offset{}, fine_y_offset{}, coarse_x_offset{}, coarse_y_offset{};
    // the scale of the coarse map compared with the smaller map.
    unsigned long scale{};
    // the length of the grid where the species start.
    long x_dim{};
    // the height of the grid where the species start.
    long y_dim{};
    double deme{};
    // for checking that the dimensions have been set before attempting to import the maps.
    bool check_set_dim;
    // for setting the movement cost through forest.
    double dispersal_relative_cost{};
    // the last time the map was updated, in generations.
    double update_time{};
    // the rate at which the habitat transforms from the modern forest map to the historical habitat map.
    // A value of 1 will give a smooth curve from the present day to historical habitat.
    double habitat_change_rate{};
    // the number of generations at which point the habitat becomes entirely historical.
    double gen_since_historical{};
    // the time the current map was updated.
    double current_map_time;
    // checks whether the simulation has already been set to the historical state.
    bool is_historical;
    // flag of whether the simulation has a historical state or not.
    bool has_historical;
    // the maximum value for habitat
    unsigned long habitat_max;
    // the maximum value on the fine map file
    unsigned long fine_max;
    // the maximum value on the coarse map file
    unsigned long coarse_max;
    // the maximum value on the historical fine map file
    unsigned long historical_fine_max;
    // the maximum value on the historical coarse map file
    unsigned long historical_coarse_max;
    // the landscape structure type
    string landscape_type;
    // true if the landscapes boundaries are infinite
    bool infinite_boundaries;
    string NextMap;
    // If this is false, there is no coarse map defined, so ignore the boundaries.
    bool has_coarse;
    // the number of updates to have occured.
    unsigned int nUpdate{};

    // Typedef for single application of the infinite landscape verses bounded landscape.
    typedef unsigned long (Landscape::*fptr)(const double &x, const double &y, const long &xwrap, const long &ywrap,
                                             const double &dCurrentGen);

    fptr getValFunc;
public:
    /**
     * @brief The default constructor.
     */
    Landscape() : mapvars(make_shared<SimParameters>())
    {
        check_set_dim = false; // sets the check to false.
        is_historical = false;
        current_map_time = 0;
        habitat_max = 1;
        getValFunc = nullptr;
        has_coarse = false;
        has_historical = false;
        landscape_type = "closed";
        fine_max = 0;
        coarse_max = 0;
        historical_fine_max = 0;
        historical_coarse_max = 0;
        infinite_boundaries = false;
    }

    /**
     * @brief Gets the maximum habitat value from any map
     * @return the maximum habitat value
     */
    unsigned long getHabitatMax();

    /**
     * @brief Returns if the simulation is using historical maps.
     * @return true if using historical maps
     */
    bool hasHistorical();

    /**
     * @brief Sets the dimensions of the grid, the area where the species are initially sampled from.
     * This function must be run before any of the calc map functions to allow for the correct deme allocation.
     *
     * @param mapvarsin the SimParameters object containing the map variables to import
     */
    void setDims(shared_ptr<SimParameters> mapvarsin);

    /**
     * @brief Checks that the map files exist (or are none/null).
     * @return true if all the paths exist in configs
     */
    bool checkMapExists();

    /**
     * @brief Imports the fine map object from file and calculates the correct values at each point.
     * Without a map to input, the fine map will simply be a matrix of 1s.
     */
    void calcFineMap();

    /**
      * @brief Imports the historical fine map object from file and calculates the correct values at each point.
      * Without a map to input, the historical fine map will simply be a matrix of 1s.
      * This has the potential to be changed easily in future versions.
      */
    void calcHistoricalFineMap();

    /**
      * @brief Imports the coarse map object from file and calculates the correct values at each point.
      * Without a map to input, the coarse map will simply be a matrix of 1s.
      * This has the potential to be changed easily in future versions.
      */
    void calcCoarseMap();

    /**
      * @brief Imports the historical coarse map object from file and calculates the correct values at each point.
      * Without a map to input, the historical coarse map will simply be a matrix of 1s.
      * This has the potential to be changed easily in future versions.
      */
    void calcHistoricalCoarseMap();

    /**
     * @brief Sets the time variables.
     * @param gen_since_historical_in the time (in generations) since a historical habitat state was achieved.
     * @param habitat_change_rate_in the rate of transform of the habitat up until the historical time.
     * A value of 0.2 would mean 20% of the change occurs linearlly up until the historical time and the remaining 80%
     * occurs in a jump to the historical state.
     */
    void setTimeVars(double gen_since_historical_in, double habitat_change_rate_in);

    /**
     * @brief Calculates the offset and extremeties of the fine map.
     *
     * Note that setting dispersal_relative_cost to a value other than 1 can massively increase simulation time.
     *
     */
    void calcOffset();

    /**
     * @brief Checks that all dimensions for all maps are zero.
     *
     * If this is true, then it means we can calculate actual offsets and dimensions from the maps, otherwise the values
     * from the parameters will be used.
     *
     * @return true if all map offsets are zero
     */
    bool checkAllDimensionsZero();

    /**
     * @brief Calculates the offsets from the map files directly.
     *
     * Assumes that all required maps have been imported.
     */
    void calculateOffsetsFromMaps();

    /**
     * @brief Uses the inputted parameters to set the offsets for the map files.
     */
    void calculateOffsetsFromParameters();

    /**
     * @brief Checks that the map file sizes are correct and that each value on the fragmented maps is less than the historical maps.
     * This should be disabled in simulations where habitat sizes are expected to shrink as well as grow.
     */
    void validateMaps();

    /**
     * @brief Checks if an update needs to be performed to the map configuration, and if it does, performs the update.
     * @param generation the current generation timer
     * @return true if the map has been updated
     */
    bool updateMap(double generation);

    /**
     * @brief Updates the historical map configuration.
     */
    void doUpdate();

    /**
     * @brief Resets the historical variables to recalculate historical maps.
     *
     * Required for rcoalescence compatability.
     */
    void resetHistorical();

    /**
     * @brief Gets the historical boolean.
     * @return the historical map state.
     */
    bool isHistorical()
    {
        if(has_historical)
        {
            return is_historical;
        }
        return true;
    }

    /**
     * @brief Sets the historical state of the system.
     * @param historical_in the historical state.
     */
    void setHistorical(const bool &historical_in)
    {
        is_historical = historical_in;
    }

    /**
     * @brief Get the historical map time
     * @return double the historical map time
     */
    double getHistorical()
    {
        return gen_since_historical;
    }

    string getLandscapeType()
    {
        return landscape_type;
    }

    /**
     * @brief Checks if the historical state has been reached.
     *
     * If there are no historical maps, this function will do nothing.
     * @param generation the time to check at.
     */
    void checkHistorical(double generation)
    {
        if(has_historical)
        {
            if(generation >= gen_since_historical)
            {
                is_historical = true;
            }
        }
    }

    /**
     * @brief Sets the landscape functions to either infinite or finite
     * @param is_infinite a string of either closed, infinite, tiled_fine or tiled_coarse,
     * corresponding to the relevant landscape type.
     *
     *
     */
    void setLandscape(string is_infinite);

    /**
     * @brief Gets the value at a particular coordinate from the correct map.
     * Takes in to account temporal and spatial referencing.
     * This version involves a call to the function pointer, *getValFunc, so that the correct call to either
     * getValFinite() or getValInfinite is made.
     * @param x the x position on the grid.
     * @param y the y position on the grid.
     * @param xwrap the number of wraps in the x dimension..
     * @param ywrap the number of wraps in the y dimension..
     * @param current_generation the current generation time.
     * @return the value on the correct map at the correct space.
     */
    unsigned long getVal(const double &x, const double &y,
                         const long &xwrap, const long &ywrap, const double &current_generation);

    /**
     * @brief Gets the value from the coarse maps, including linear interpolating between the historical and present maps
     * @param xval the x coordinate
     * @param yval the y coordinate
     * @param current_generation the current generation timer
     * @return the value of the map at the given coordinates and time
     */
    unsigned long getValCoarse(const double &xval, const double &yval, const double &current_generation);

    /**
     * @brief Gets the value from the fine maps, including linear interpolating between the historical and present maps
     * @param xval the x coordinate
     * @param yval the y coordinate
     * @param current_generation the current generation timer
     * @return the value of the map at the given coordinates and time
     */
    unsigned long getValFine(const double &xval, const double &yval, const double &current_generation);

    /**
     * @brief Gets the value at a particular coordinate from the correct map.
     * Takes in to account temporal and spatial referencing. This version assumes finite landscape.
     * @param x the x position on the grid.
     * @param y the y position on the grid.
     * @param xwrap the number of wraps in the x dimension..
     * @param ywrap the number of wraps in the y dimension..
     * @param current_generation the current generation time.
     * @return the value on the correct map at the correct space.
     */
    unsigned long getValFinite(const double &x, const double &y, const long &xwrap, const long &ywrap,
                               const double &current_generation);

    /**
     * @brief Gets the value at a particular coordinate from the correct map.
     * Takes in to account temporal and spatial referencing. This version assumes an infinite landscape.
     * @param x the x position on the grid.
     * @param y the y position on the grid.
     * @param xwrap the number of wraps in the x dimension..
     * @param ywrap the number of wraps in the y dimension..
     * @param current_generation the current generation time.
     * @return the value on the correct map at the correct space.
     */
    unsigned long getValInfinite(const double &x, const double &y, const long &xwrap, const long &ywrap,
                                 const double &current_generation);

    /**
     * @brief Gets the value at a particular coordinate from the correct map.
     * Takes in to account temporal and spatial referencing.
     * This version assumes an infinite landscape of tiled coarse maps.
     * @param x the x position on the grid.
     * @param y the y position on the grid.
     * @param xwrap the number of wraps in the x dimension..
     * @param ywrap the number of wraps in the y dimension..
     * @param current_generation the current generation time.
     * @return the value on the correct map at the correct space.
     */
    unsigned long getValCoarseTiled(const double &x, const double &y, const long &xwrap, const long &ywrap,
                                    const double &current_generation);

    /**
     * @brief Gets the value at a particular coordinate from the correct map.
     * Takes in to account temporal and spatial referencing.
     * This version assumes an infinite landscape of tiled fine maps.
     * @param x the x position on the grid.
     * @param y the y position on the grid.
     * @param xwrap the number of wraps in the x dimension..
     * @param ywrap the number of wraps in the y dimension..
     * @param current_generation the current generation time.
     * @return the value on the correct map at the correct space.
     */
    unsigned long getValFineTiled(const double &x, const double &y, const long &xwrap, const long &ywrap,
                                  const double &current_generation);

    /**
     * @brief Gets the x position on the fine map, given an x and x wrapping.
     *
     * Note that this function will not check if the value is actually within bounds of the fine map,
     * and an error will likely be thrown by the matrix referencing if this is the case.
     * @param x the x coordinate on the sample mask
     * @param xwrap the x wrapping of the sample mask.
     * @return the x location on the fine map
     */
    unsigned long convertSampleXToFineX(const unsigned long &x, const long &xwrap);

    /**
     * @brief Gets the y position on the fine map, given a y and y wrapping.
     *
     * Note that this function will not check if the value is actually within bounds of the fine map,
     * and an error will likely be thrown by the matrix referencing if this is the case.
     * @param y the y coordinate on the sample mask
     * @param ywrap the y wrapping of the sample mask.
     * @return the y location on the fine map
     */
    unsigned long convertSampleYToFineY(const unsigned long &y, const long &ywrap);

    /**
     * @brief Converts the fine map coordinates to the sample grid coordinates.
     * Main conversion is in a call to convertCoordinates, but also makes sure the returned types are long integers.
     * @param x the x coordinate to modify
     * @param xwrap the x wrapping to modify
     * @param y the y coordinate to modify
     * @param ywrap the y wrapping to modify
     */
    void convertFineToSample(long &x, long &xwrap, long &y, long &ywrap);

    /**
     * @brief Counts the number of spaces available in the initial species space. Requires the samplemask to check the sampling area.
     * @param dSample the sample proportion (from 0 to 1).
     * @param samplemask the DataMask object to sample from.
     * @return the total number of individuals predicted to initially exist on the map.
     */
    unsigned long getInitialCount(double dSample, DataMask &samplemask);

    /**
     * @brief Gets the mapvars object pointer for referencing simulation parameters.
     * @return
     */
    shared_ptr<SimParameters> getSimParameters();

    /**
     * @brief Checks whether the point is habitat or non-habitat.
     *@param x the x position on the grid.
     * @param y the y position on the grid.
     * @param xwrap the number of wraps in the x dimension.
     * @param ywrap the number of wraps in the y dimension.
     * @param generation the current generation time.
     * @return a boolean of whether the map is habitat or non-habitat.
     */
    bool checkMap(const double &x, const double &y, const long &xwrap, const long &ywrap, const double &generation);

    /**
     * @brief  Checks whether the point comes from the fine grid.
     * @param x the x position
     * @param y the y position
     * @param xwrap the number of wraps in the x dimension
     * @param ywrap the number of wraps in the y dimension
     * @return a boolean of whether the location is on the fine map
     */
    bool isOnFine(const double &x, const double &y, const long &xwrap, const long &ywrap);

    /**
     * @brief  Checks whether the point comes from the coarse grid.
     * @param x the x position
     * @param y the y position
     * @param xwrap the number of wraps in the x dimension
     * @param ywrap the number of wraps in the y dimension
     * @return a boolean of whether the location is on the fine map
     */
    bool isOnCoarse(const double &x, const double &y, const long &xwrap, const long &ywrap);

    /**
     * @brief Checks that the point supplied is within map limits.
     * If the map is inifite, returns true.
     * @param x the x position
     * @param y the y position
     * @param xwrap the number of wraps in the x dimension
     * @param ywrap the number of wraps in the y dimension
     * @return a boolean of whether the location is on the fine map
     * @return true if the point is within the map limits
     */
    bool isOnMap(const double &x, const double &y, const long &xwrap, const long &ywrap);

    /**
     * @brief Fixes the coordinates to be correctly within the original grid, altering the xwrap and ywrap consequently.
     * @param x the x position.
     * @param y the y position.
     * @param xwrap the number of wraps in the x dimension.
     * @param ywrap the number of wraps in the y dimension.
     */
    void fixGridCoordinates(double &x, double &y, long &xwrap, long &ywrap);

    /**
     * @brief The function that actually performs the dispersal.
     * It is included here for easier  programming and efficiency as the function doesn't need to perform all the checks
     * until the edge of the fine grid.
     * @param dist the distance travelled (or "distance energy" if dispersal_relative_cost is not 1).
     * @param angle the angle of movement.
     * @param startx the start x position.
     * @param starty the start y position.
     * @param startxwrap the start number of wraps in the x dimension.
     * @param startywrap the start number of wraps in the y dimension.
     * @param disp_comp a boolean of whether the dispersal was complete or not. This value is returned false if dispersal
     * is to habitat, false otherwise.
     * @param generation the time in generations since the start of the simulation.
     * @return the density value at the end dispersal point
     */
    unsigned long runDispersal(const double &dist, const double &angle, long &startx, long &starty, long &startxwrap,
                               long &startywrap, bool &disp_comp, const double &generation);

    /**
     * @brief Calculates the distance from the start position to the nearest habitat cell.
     * @param start_x the start x coordinate
     * @param start_y the start y coordinate
     * @param start_x_wrap the starting x wrapping
     * @param start_y_wrap the starting y wrapping
     * @param generation the generation timer
     * @return the distance from the start position to the nearest habitat cell
     */
    double distanceToNearestHabitat(const long &start_x, const long &start_y, const long &start_x_wrap,
                                    const long &start_y_wrap, const double &generation);

    /**
     * @brief Gets the nearest habitat cells from a particular point, spiraling outwards.
     * @param start_x the start x coordinate
     * @param start_y the start y coordinate
     * @param start_x_wrap the starting x wrapping
     * @param start_y_wrap the starting y wrapping
     * @param end_x the end x coordinate value to modify
     * @param end_y the end y coordinate value to modify
     * @param generation the generation timer
     */
    void findNearestHabitatCell(const long &start_x, const long &start_y, const long &start_x_wrap,
                                const long &start_y_wrap, double &end_x, double &end_y, const double &generation);

    /**
     * @brief Finds the nearest habitat cell using a much slower method (scanning the entire map for cells.
     * @param start_x the start x coordinate
     * @param start_y the start y coordinate
     * @param start_x_wrap the starting x wrapping
     * @param start_y_wrap the starting y wrapping
     * @param end_x the end x coordinate value to modify
     * @param end_y the end y coordinate value to modify
     * @param generation the generation timer
     * @return true if a habitat cell is found
     */
    bool findAnyHabitatCell(const long &start_x, const long &start_y, const long &start_x_wrap,
                            const long &start_y_wrap, double &end_x, double &end_y, const double &generation);

    /**
     * @brief Operator for outputting the Map object variables to an output stream.
     * This is used for storing the Map object to file.
     * @param os the output stream.
     * @param r the Map object to output.
     * @return the output stream.
     */
    friend ostream &operator<<(ostream &os, const Landscape &r)
    {
        os << r.fine_x_min << "\n" << r.fine_x_max << "\n" << r.coarse_x_min << "\n"
           << r.coarse_x_max;
        os << "\n" << r.fine_y_min << "\n" << r.fine_y_max << "\n" << r.coarse_y_min << "\n" << r.coarse_y_max << "\n";
        os << r.fine_x_offset << "\n" << r.fine_y_offset << "\n" << r.coarse_x_offset << "\n" << r.coarse_y_offset
           << "\n";
        os << r.scale << "\n" << r.x_dim << "\n" << r.y_dim << "\n" << r.deme << "\n" << r.check_set_dim << "\n"
           << r.dispersal_relative_cost << "\n";
        os << r.update_time << "\n" << r.habitat_change_rate << "\n" << r.gen_since_historical << "\n"
           << r.current_map_time << "\n"
           << r.is_historical << "\n";
        os << r.NextMap << "\n" << r.nUpdate << "\n" << r.landscape_type << "\n" << r.fine_max << "\n"
           << r.coarse_max << "\n";
        os << r.historical_fine_max << "\n" << r.historical_coarse_max << "\n" << r.habitat_max << "\n"
           << r.has_coarse << "\n" << r.has_historical << "\n";
        return os;
    }

    /**
     * @brief Operator for inputting the Map object variables from an input stream.
     * This is used for reading the Map object from file.
     * @param is the input stream.
     * @param r the Map object to input to.
     * @return the input stream.
     */
    friend istream &operator>>(istream &is, Landscape &r)
    {
        is >> r.fine_x_min;
        is >> r.fine_x_max >> r.coarse_x_min;
        is >> r.coarse_x_max >> r.fine_y_min >> r.fine_y_max;
        is >> r.coarse_y_min >> r.coarse_y_max;
        is >> r.fine_x_offset >> r.fine_y_offset >> r.coarse_x_offset >> r.coarse_y_offset >> r.scale >> r.x_dim
           >> r.y_dim
           >> r.deme >> r.check_set_dim >> r.dispersal_relative_cost;
        is >> r.update_time >> r.habitat_change_rate >> r.gen_since_historical >> r.current_map_time >> r.is_historical;
        getline(is, r.NextMap);
        is >> r.nUpdate;
        is >> r.landscape_type;
        is >> r.fine_max >> r.coarse_max;
        is >> r.historical_fine_max >> r.historical_coarse_max;
        is >> r.habitat_max >> r.has_coarse >> r.has_historical;
        r.setLandscape(r.mapvars->landscape_type);
        r.calcFineMap();
        r.calcCoarseMap();
        r.calcHistoricalFineMap();
        r.calcHistoricalCoarseMap();
        r.recalculateHabitatMax();
        return is;
    }

    /**
     * @brief Prints some selected Map variables to the terminal.
     * @return the string containing the map variables to print
     */
    string printVars();

    /**
     * @brief Wipes the map of all variables. Only really useful for testing purposes.
     */
    void clearMap();

    /**
     * @brief Recalculates the habitat map maximum by checking the maximums for each of the relevant map files
     * (fine, coarse and historicals).
     */
    void recalculateHabitatMax();

};

#endif // LANDSCAPE_H
 