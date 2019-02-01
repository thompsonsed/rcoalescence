//This file is part of necsim project which is released under MIT license.
//See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.
/**
 * @author Samuel Thompson
 * @date 30/08/2016
 * @file Datapoint.h
 *
 * @brief Contains the DataPoint class for storing objects during simulation run time.
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 * This class is only used during simulation runs and is not outputted to a database.
 * A Row of DataPoint objects is utilised by the main Tree objects.
 */

#ifndef DATAPOINT_H
#define DATAPOINT_H

#include <iostream>
#include "Logging.h"

using namespace std;

/**
 * @brief A data object used in coalescence simulations for calculating the output.
 * Data from this object is outputted to an SQLite database after simulations are complete.
 */
class DataPoint
{

private:
    // x position
    unsigned long xpos;
    // y position
    unsigned long ypos;
    // number of wraps of x around the torus
    long xwrap;
    // number of wraps of y around the torus
    long ywrap;
    // the next individual in the loop of those that have the same xypos
    unsigned long next_lineage;
    // points to the position in output of this lineage
    unsigned long reference;
    // points to the position in the SpeciesList file.
    unsigned long list_position;
    // the reference number within the linked species_id_list of wrapped lineages
    unsigned long nwrap;
    // the max-min number
    double min_max;
public:

    /**
     * @brief Standard constructor
     */
    DataPoint() : xpos(0), ypos(0), xwrap(0), ywrap(0), next_lineage(0), reference(0), list_position(0), nwrap(0),
                  min_max(0)
    {

    }

    /**
     * @brief Standard destructor
     */
    ~DataPoint() = default;

    /**
     * @brief Setup of lineage data with any information that's wanted.
     * Note that nwrap is set to 0 in this routine.
     * @param x the x position on the grid
     * @param y the y position on the grid
     * @param xwrap_in the number of wraps of the location on the grid in the x direction
     * @param ywrap_in the number of wraps of the location on the grid in the y direction
     * @param reference_in the position in the TreeNode reference object
     * @param list_position_in the position within the SpeciesList object at the relevant x,y position
     * @param min_max_in the input maximum minimum speciation rate required for speciation to have occured
     */
    void setup(unsigned long x, unsigned long y, long xwrap_in, long ywrap_in, unsigned long reference_in,
               unsigned long list_position_in, double min_max_in);

    /**
     * @brief Setup of lineage data with any information that's wanted.
     * @param reference_in the reference point for the TreeNode reference object
     * @param list_position_in the list position of this lineage in SpeciesList object
     * @param min_max_in the input maximum minimum speciation rate required for speciation to have occured
     */
    void setup(unsigned long reference_in, unsigned long list_position_in, double min_max_in);

    /**
     * @brief Copy constructor from another Datapoint object.
     * @param datin a Datapoint object to copy the data from.
     */
    void setup(DataPoint datin);

    /**
     * @brief Sets the mpos (the position within the Row of TreeNode objects.).
     * @param z the desired mpos.
     */
    void setReference(unsigned long z);

    /**
     * @brief Set the next link in the linked list.
     * @param x the next Datapoint object.
     */
    void setNext(unsigned long x);

    /**
     * @brief Sets the list position within the SpeciesList object.
     * @param l the input list position.
     */
    void setListPosition(unsigned long l);

    /**
     * @brief Sets the number of wraps from the first SpeciesList wrapped lineage.
     * If and only if this is 0, the lineage is within the main grid (i.e xwrap and ywrap should be 0).
     * @param n the desired nwrap.
     */
    void setNwrap(unsigned long n);

    /**
     * @brief Sets the minmax variable.
     * This is the minimum maximum speciation rate required for speciation to have occured on this branch.
     * @param d the minmax to set.
     */
    void setMinmax(double d);

    /**
     * @brief Get the x position.
     * @return the xpos.
     */
    unsigned long getXpos();

    /**
     * @brief Get the y position.
     * @return the ypos.
     */
    unsigned long getYpos();

    /**
     * @brief Get the x wrapping.
     * @return the xwrap.
     */
    long getXwrap();

    /**
     * @brief Get the y wrapping.
     * @return the ywrap.
     */
    long getYwrap();

    /**
     * @brief Get the reference position variable.
     * @return the mpos.
     */
    unsigned long getReference();

    /**
     * @brief Gets the next element linked to this DataPoint
     * @return the reference of the next individual in the linked list
     */
    unsigned long getNext();

    /**
     * @brief Gets the list position with the SpeciesList object at the relevant x,y position.
     * @return the listpos.
     */
    unsigned long getListpos();

    /**
     * @brief Get the position in the linked list from the SpeciesList object.
     * If this is 0, indicates the lineage lies on the original grid, and xwrap and ywrap should be 0.
     * @return the nwrap.
     */
    unsigned long getNwrap();

    /**
     * @brief Get the maximum minimum speciation rate required for speciation to have occured on this branch.
     * @return the minmax.
     */
    double getMinmax();

    /**
     * @brief Decreases the nwrap by 1 (to a minimum of 0).
     */
    void decreaseNwrap();

    /**
     * @brief Sets the position in space.
     * @param x the x position.
     * @param y the y position.
     * @param xwrapin the number of wraps in the x direction.
     * @param ywrapin the number of wraps in the y direction.
     */
    void setEndpoint(long x, long y, long xwrapin, long ywrapin);

    /**
     * @brief An operator for piping the variables of the Datapoint object to the output stream.
     * @param os the output stream.
     * @param d the Datapoint object to output.
     * @return returns the output stream at the end.
     */
    friend ostream &operator<<(ostream &os, const DataPoint &d);

    /**
     * @brief An operator for piping the variables in to the Datapoint object from the input stream.
     * @param is the input stream
     * @param d the Datapoint object to input to.
     * @return returns the input stream at the end.
     */
    friend istream &operator>>(istream &is, DataPoint &d);

#ifdef DEBUG
    void logActive(const int &level);
#endif // DEBUG
};

#endif // DATAPOINT_H