// This file is part of necsim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @author Samuel Thompson
 * @file SpeciesAbundancesHandler.h
 *
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 * @brief Base class for storing and generating species abundances
 *
 * Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
 */

#ifndef SPECIES_ABUNDANCES_H
#define SPECIES_ABUNDANCES_H

#include <vector>
#include <map>
#include <memory>
#include "RNGController.h"

using namespace std;

/**
 * @brief Base class for defining species abundances.
 */
class SpeciesAbundancesHandler
{
protected:

    shared_ptr<RNGController> random;
    unsigned long max_species_id;
    unsigned long metacommunity_size;
    unsigned long local_community_size;
    long double speciation_rate;
public:

    /**
     * @brief Default constructor
     */
    SpeciesAbundancesHandler();

    /**
     * @brief Default destructor
     */
    virtual ~SpeciesAbundancesHandler() = default;;

    /**
     * @brief Creates the SpeciesAbundancesHandler object
     * @param random the random number generator
     * @param metacommunity_size the number of individuals in the metacommunity
     * @param speciation_rate the speciation rate of the metacommunity
     * @param local_community_size: the number of individuals in the local community
     */
    virtual void setup(shared_ptr<RNGController> random, const unsigned long &metacommunity_size,
                       const long double &speciation_rate,
                       const unsigned long &local_community_size);

    /**
     * @brief Gets a randomly generated species identity.
     * @return the species identity
     */
    virtual unsigned long getRandomSpeciesID() = 0;

    /**
     * @brief Sets the abundance list.
     * @param abundance_list_in list of abundances for each species
     */
    virtual void setAbundanceList(const shared_ptr<map<unsigned long, unsigned long>> &abundance_list_in);

    /**
     * @brief Sets the abundance list.
     * @param abundance_list_in list of abundances for each species
     */
    virtual void setAbundanceList(shared_ptr<vector<unsigned long>> abundance_list_in);
};

#endif //SPECIES_ABUNDANCES_H
