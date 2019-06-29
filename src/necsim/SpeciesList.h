//This file is part of necsim project which is released under MIT license.
//See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @author Samuel Thompson
 * @file SpeciesList.h
 * @brief  Contains the SpeciesList class for usage in coalescence simulations.
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */

#ifndef SPECIESLIST
#define SPECIESLIST

/**
 * @class SpeciesList
 * @brief Contains a list of the species that exist at one location.
 * The Row object, list, contains the active reference number, for looking up the lineage in a Row of Datapoint objects.
 * Also contains the functions for correctly generating coalescence probabilities and list management.
 * 
 * Note that the maximum size of the list is constrained by the maximum size of unsigned long. Any simulation requiring
 * more individuals per cell than this will unlikely finish in any reasonable time anyway.
 */
#ifdef CXX14_SUPPORT
#include "memory.h"
#else

#include <memory>

#endif

#include "Matrix.h"
#include "RNGController.h"

using namespace std;

class SpeciesList
{
private:
    unsigned long list_size, maxsize; // List size and maximum size of the cell (based on percentage cover).
    unsigned long next_active; // For calculating the wrapping, using the next and last system.
    vector<unsigned long> species_id_list; // list of the active reference number, with zeros for empty cells.
    unsigned long nwrap; // The number of wrapping (next and last possibilities) that there are.
public:
    /**
     * @brief Default constructor
     */
    SpeciesList();

    /**
     * @brief Default destructor
     */
    ~SpeciesList() = default;

    /**
     * @brief Initialises the list to the specified size.
     * @param maxsizein the maximum list size.
     */
    void initialise(unsigned long maxsizein);

    // special case if just the maxsize wants to be change, but want to maintain the species_id_list variables.
    /**
     * @brief Sets the maxsize without altering the actual size of list.
     * @param maxsizein The new maximum size to set.
     */
    void setMaxsize(unsigned long maxsizein);

    /**
     * @brief Set specific entry to a particular species reference number.
     * @param index the location in list of the species.
     * @param new_val the new species reference to set list[index] to.
     *
     * This version throws an error if the space is empty.
     */
    void setSpecies(unsigned long index, unsigned long new_val);

    /**
     * @brief Set specific entry to a particular species reference number
     * @param index the location in list of the species
     * @param new_val the new species reference to set list[index] to
     *
     * Note this version will throw a runtime_error if the space is not empty
     */
    void setSpeciesEmpty(unsigned long index, unsigned long new_val);

    /**
     * @brief Set the next active lineage (for wrapping purposes).
     * @param n the lineage to set as the first wrapped lineage.
     */
    void setNext(unsigned long n);

    /**
     * @brief Set the number of wrapping lineages.
     * @param nr the number of wrapped lineages.
     */
    void setNwrap(unsigned long nr);

    /**
     * @brief Add a new species to the first empty place and return the position of the lineage.
     * @param new_spec the new species reference to place in the first empty space.
     * @return the location the species has been added to.
     */
    unsigned long addSpecies(const unsigned long &new_spec);

    /**
     * @brief Removes the species at the specified index.
     * The species number will be replaced with 0, indicating no species present.
     *
     * Older versions of this function re-shuffled the list so that all species came at the top.
     * @param index the index of the species to remove from the list.
     */
    void deleteSpecies(unsigned long index);

    /**
     * @brief Decreases the nwrap by one.
     *
     * Indicates the number of species wrapped at the location of this SpeciesList object has decreased by one.
     */
    void decreaseNwrap();

    /**
     * @brief Increases the list size by one.
     */
    void increaseListSize();

    /**
     * @brief Increases the nwrap by one.
     *
     * Indicates the number of species wrapped at the location of this SpeciesList object has increased by one.
     */
    void increaseNwrap();

    /**
     * @brief Changes the maximum size of the SpeciesList.
     * Currently identical to setMaxsize.
     * @param newmaxsize the new maximum size to be applied.
     */
    void changePercentCover(unsigned long newmaxsize);

    /**
     * @brief Get a random species reference number from all the potential entries.
     * Updated alternative version returns any entry, including empty cells, giving the probability of coalescence as well.
     * @param rand_no the random number object to pass (for maintaining the same seed throughout simulations).
     * @return the reference of the random lineage. 0 indicates an empty space.
     */
    unsigned long getRandLineage(shared_ptr<RNGController> rand_no);

    /**
     * @brief Get the species reference number from a particular entry.
     * @param index the location of the species to reference.
     * @return the species reference at the specified location.
     */
    unsigned long getSpecies(unsigned long index);

    /**
     * @brief Get the next_active variable.
     * @return the next linked species reference.
     */
    unsigned long getNext();

    /**
     * @brief Getter for the nwrap.
     * @return the number of wrapped lineages currently at this grid cell.
     */
    unsigned long getNwrap();

    /**
     * @brief Getter for the list size.
     * @return the number of lineages currently directly within the SpeciesList.
     */
    unsigned long getListSize();

    /**
     * @brief Getter for the maximum size of the SpeciesList object.
     * @return the maximum number of lineages that can exist currently.
     */
    unsigned long getMaxSize();

    /**
     * @brief Gets the length of the list object
     * @return the length of the list object
     */
    unsigned long getListLength();

    /**
     * @brief Empties the list of any data and fills the list with zeros.
     */
    void wipeList();

    /**
     * @brief Outputs the SpeciesList object to an output stream.
     * Allows for piping to the terminal or writing the object to a file.
     * @param os the output stream.
     * @param r the SpeciesList object to output.
     * @return the output stream.
     */
    friend ostream &operator<<(ostream &os, const SpeciesList &r);

    /**
     * @brief Inputs the SpeciesList object from an input stream.
     * Allows for reading data from a file or string stream.
     * @param is the input stream.
     * @param r the SpeciesList object to input to.
     */
    friend istream &operator>>(istream &is, SpeciesList &r);
};

#endif