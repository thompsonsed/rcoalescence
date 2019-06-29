//This file is part of necsim project which is released under MIT license.
//See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @author Samuel Thompson
 * @date 07/08/2017
 * @file DispersalCoordinator.h
 * @brief Contains the DispersalCoordinator, which contains all routines related to dispersal
 * including utilisation of density maps and dispersal probability maps.
 * 
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */
#ifndef DISPERSALCOORDINATOR_H
#define DISPERSALCOORDINATOR_H

#include <cstring>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <cmath>
#include <memory>

#include "RNGController.h"
#include "Map.h"
#include "Step.h"
#include "Landscape.h"
#include "ActivityMap.h"

/**
 * @brief Class for generating dispersal distances and provide routines for reading dispersal distance maps
 * as a unwound map-of-maps. This class also handles reading density maps for rejection sampling.
 * 
 * It requires linking to a density map, random number generator and a generation counter from the Tree class.
 * 
 * Note that no element of this object is recorded during a paused simulation, as all objects pointed to are stored 
 * elsewhere and behaviours are recalculated upon simulation resume.
 */
class DispersalCoordinator
{
protected:

    // Our map of dispersal probabilities (if required)
    // This will contain cummulative probabilities across rows
    // So dispersal is from the y cell to each of the x cells.
    Map<double> dispersal_prob_map;
    // This object is only used if there are multiple density maps over time.
    Map<double> raw_dispersal_prob_map;
    // Our random number generator for dispersal distances
    // This is a pointer so that the random number generator is the same
    // across the program.

    shared_ptr<RNGController> NR;
    // Pointer to the landscape object for getting density values.
    shared_ptr<Landscape> landscape;
    // Pointer to the reproduction map object for obtaining reproduction probabilities
    shared_ptr<ActivityMap> reproduction_map;
    // Pointer to the generation counter for the simulation
    double *generation;

    // function ptr for our getDispersal function
    typedef void (DispersalCoordinator::*dispersal_fptr)(Step &this_step);

    dispersal_fptr doDispersal;

    // Function pointer for end checks
    typedef bool (DispersalCoordinator::*end_fptr)(const unsigned long &density, long &oldx, long &oldy, long &oldxwrap,
                                                   long &oldywrap, const long &startx, const long &starty,
                                                   const long &startxwrap, const long &startywrap);

    // once setup will contain the end check function to use for this simulation.
    end_fptr checkEndPointFptr;
    unsigned long xdim;
    unsigned long ydim;

public:
    DispersalCoordinator();

    ~DispersalCoordinator();

    /**
     * @brief Sets the random number pointer to an NRrand instance.
     * @param NR_ptr the random number object to set to
     */
    void setRandomNumber(shared_ptr<RNGController> NR_ptr);

    /**
     * @brief Sets the pointer to the Landscape object
     * @param landscape_ptr pointer to a Landscape object
     * @param repr_map_ptr pointer to the reproduction probability map
     */
    void setMaps(shared_ptr<Landscape> landscape_ptr, shared_ptr<ActivityMap> repr_map_ptr);

    /**
     * @brief Sets the pointer to the Landscape object
     * @param landscape_ptr pointer to a Landscape object
     */
    void setMaps(shared_ptr<Landscape> landscape_ptr);

    /**
     * @brief Sets the generation pointer to the provided double
     * @param generation_ptr pointer to the generation double
     */
    void setGenerationPtr(double *generation_ptr);

    /**
     * @brief Sets the dispersal method and parameters
     *
     * @param dispersal_method string containing the dispersal type. Can be one of [normal, fat-tail, norm-uniform]
     * @param dispersal_file string containing the dispersal file, or "none" if using dispersal kernel
     * @param dispersal_x the x dimensions of the dispersal file
     * @param dispersal_y the y dimensions of the dispersal file
     * @param m_probin the probability of drawing from the uniform distribution. Only relevant for uniform dispersals
     * @param cutoffin the maximum value to be drawn from the uniform dispersal. Only relevant for uniform dispersals
     * @param sigmain the fatness of the fat-tailed dispersal kernel
     * @param tauin the width of the fat-tailed dispersal kernel
     * @param restrict_self if true, denies possibility that dispersal comes from the same cell as the parent
     */
    void setDispersal(const string &dispersal_method, const string &dispersal_file,
                      const unsigned long &dispersal_x, const unsigned long &dispersal_y,
                      const double &m_probin, const double &cutoffin,
                      const double &sigmain, const double &tauin, const bool &restrict_self);

    /**
     * @brief Sets the dispersal parameters from the SimParameters object.
     * @param simParameters pointer to the simulation parameters to set
     */
    void setDispersal(shared_ptr<SimParameters> simParameters);

    /**
     * @brief Imports the dispersal map and fixes the values based on the density and reproduction probabilities.
     * @param dispersal_dim pointer x and y dimension of the dispersal map
     * @param dispersal_file name of the dispersal file
     */
    void importDispersal(const unsigned long &dispersal_dim, const string &dispersal_file);

    /**
     * @brief Saves the raw dispersal data to another memory object to save reading from disk multiple times later on.
     * This is slightly more RAM intensive which may cause issues for certain simulations.
     */
    void setRawDispersalMap();

    /**
     * @brief Adds the density values to the dispersal map.
     * @param generation the current generation counter
     */
    void addDensity();

    /**
     * @brief Adds the reproduction rates to the dispersal map.
     */
    void addReproduction();

    /**
     * @brief Fixes the dispersal map by generating cumulative probability distributions across each row.
     */
    void fixDispersal();

    /**
     * @brief Sums probabilities across the given row and divides by the total probability, so that the cumulative
     * probabilities sum to one.
     * @param row the row index to check over
     */
    void fixDispersalRow(unsigned long row);

    /**
     * @brief Checks if the provided row index requires re-scaling to a cumulative probability in the dispersal map.
     * @param row the row index to check
     * @return true if re-scaling of the row is required
     */
    bool checkDispersalRow(unsigned long row);

    /**
     * @brief Ensures that the dispersal map makes sense given the density.
     */
    void verifyDispersalMapDensity();

    /**
     * @brief Checks that the dispersal map makes sense with dispersal to and from only cells which have a non-zero
     * density.
     */
    void verifyDispersalMapSetup();

    /**
     * @brief Updates the dispersal map, if there is one, to reflect changes in landscape density.
     */
    void updateDispersalMap();

#ifdef DEBUG
    /**
     * @brief Asserts that the cell reference is correct for the provided coordinates.
     * @param expected the expected cell reference.
     */
    void assertReferenceMatches(unsigned long expected);
#endif // DEBUG

    /**
     * @brief Picks a random cell from the whole map and stores the value in the step object
     * @param this_step the step object to store end points in
     */
    void disperseNullDispersalMap(Step &this_step);

    /**
     * @brief Picks a random dispersal distance from the dispersal map
     * @param this_step the step object to store end points in
     */
    void disperseDispersalMap(Step &this_step);

    /**
     * @brief Calculates the new coordinates for a column reference.
     * This includes converting between the fine map and sample map.
     * New coordinates are saved in this_step.
     * @param this_step the step to save new coordinates in.
     * @param col_ref the column reference for
     */
    void calculateCellCoordinates(Step &this_step, const unsigned long &col_ref);

    /**
     * @brief Calculates the cell reference for a particular coordinate
     *
     * The formula for this calculation is x + (y * xdim) where xdim is the x dimension of the fine map,
     * and x and y are the coordinates for the fine map
     *
     * @param this_step the step object containing the x, y location, and x,y wrapping
     * @return the cell reference from the dispersal_prob_map which corresponds to the required cell
     */
    unsigned long calculateCellReference(Step &this_step);

    /**
     * @brief Calls the dispersal kernel from the supplied dispersal distribution.
     * @param this_step the step object to store end points in
     */
    void disperseDensityMap(Step &this_step);

    /**
     * @brief Disperses to the nearest habitat cell from a given location.
     * @param this_step teh step object to store end points in
     */
    void disperseNearestHabitat(Step &this_step);

    /**
     * @brief Sets the end point function pointer correctly, based on whether it is restricted or not.
     * @param restrict_self if true, denies possibility that dispersal comes from the same cell as the parent
     */
    void setEndPointFptr(const bool &restrict_self);

    /**
     * @brief Check the end point for the given coordinates and density
     * @param density the density at the end point - avoids an extra call to Map::getVal()
     * @param oldx the old x position
     * @param oldy the old y position
     * @param oldxwrap the old x wrap
     * @param oldywrap the old y wrap
     * @param startx the starting x position
     * @param starty the starting y position
     * @param startxwrap the starting x wrap
     * @param startywrap the ending y wrap
     * @return true if the end point passes the density and restricted checks
     */
    bool checkEndPoint(const unsigned long &density, long &oldx, long &oldy, long &oldxwrap, long &oldywrap,
                       const long &startx, const long &starty, const long &startxwrap, const long &startywrap);

    /**
     * @brief Check the end point for the given coordinates and density
     * @param density the density at the end point - avoids an extra call to Map::getVal()
     * @param oldx the old x position
     * @param oldy the old y position
     * @param oldxwrap the old x wrap
     * @param oldywrap the old y wrap
     * @param startx the starting x position
     * @param starty the starting y position
     * @param startxwrap the starting x wrap
     * @param startywrap the ending y wrap
     * @return true if the end point passes the density and restricted checks
     */
    bool checkEndPointDensity(const unsigned long &density, long &oldx, long &oldy, long &oldxwrap, long &oldywrap,
                              const long &startx, const long &starty, const long &startxwrap, const long &startywrap);

    /**
     * @brief Check the end point for the given coordinates and density
     * @param density the density at the end point - avoids an extra call to Map::getVal()
     * @param oldx the old x position
     * @param oldy the old y position
     * @param oldxwrap the old x wrap
     * @param oldywrap the old y wrap
     * @param startx the starting x position
     * @param starty the starting y position
     * @param startxwrap the starting x wrap
     * @param startywrap the ending y wrap
     * @return true if the end point passes the density and restricted checks
     */
    bool checkEndPointRestricted(const unsigned long &density, long &oldx, long &oldy, long &oldxwrap, long &oldywrap,
                                 const long &startx, const long &starty, const long &startxwrap,
                                 const long &startywrap);

    /**
     * @brief Check the end point for the given coordinates and density
     * @param density the density at the end point - avoids an extra call to Map::getVal()
     * @param oldx the old x position
     * @param oldy the old y position
     * @param oldxwrap the old x wrap
     * @param oldywrap the old y wrap
     * @param startx the starting x position
     * @param starty the starting y position
     * @param startxwrap the starting x wrap
     * @param startywrap the ending y wrap
     * @return true if the end point passes the density and restricted checks
     */
    bool checkEndPointDensityReproduction(const unsigned long &density, long &oldx, long &oldy, long &oldxwrap,
                                          long &oldywrap, const long &startx, const long &starty,
                                          const long &startxwrap, const long &startywrap);

    /**
     * @brief Check the end point for the given coordinates and density
     * @param density the density at the end point - avoids an extra call to Map::getVal()
     * @param oldx the old x position
     * @param oldy the old y position
     * @param oldxwrap the old x wrap
     * @param oldywrap the old y wrap
     * @param startx the starting x position
     * @param starty the starting y position
     * @param startxwrap the starting x wrap
     * @param startywrap the ending y wrap
     * @return true if the end point passes the density and restricted checks
     */
    bool checkEndPointDensityRestrictedReproduction(const unsigned long &density, long &oldx, long &oldy,
                                                    long &oldxwrap, long &oldywrap, const long &startx,
                                                    const long &starty, const long &startxwrap, const long &startywrap);

    /**
     * @brief Performs the dispersal routine using the Step object to read starting positions and record the end positions.
     * @param this_step the Step object for reading starting position and storing output distances and angles
     */
    void disperse(Step &this_step);

};

#endif // DISPERSAL_H
