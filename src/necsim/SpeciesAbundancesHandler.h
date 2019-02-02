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
#include "NRrand.h"

using namespace std;

/**
 * @brief Base class for defining species abundances.
 */
class SpeciesAbundancesHandler
{
protected:

    shared_ptr<NRrand> random;
    unsigned long max_species_id;
    unsigned long community_size;
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
     * @param community_size the number of individuals in the community
     * @param speciation_rate the speciation rate of the community
     */
    virtual void setup(shared_ptr<NRrand> random, const unsigned long &community_size,
                       const long double &speciation_rate);

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

    /**
     * @brief Gets a random species abundance.
     * @return the randomly generated abundance
     */
    virtual unsigned long getRandomAbundanceOfIndividual();

    /**
     * @brief Gets the species richness of a particular abundance class.
     * @param abundance the abundance class of the species
     * @return the number of species with that abundance
     */
    virtual unsigned long getSpeciesRichnessOfAbundance(const unsigned long &abundance)
    {
        return 0;
    }

};

#endif //SPECIES_ABUNDANCES_H
