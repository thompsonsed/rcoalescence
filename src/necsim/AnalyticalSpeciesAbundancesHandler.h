// This file is part of necsim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @author Samuel Thompson
 * @file AnalyticalSpeciesAbundancesHandler.h
 *
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 * @brief Class for repeatedly selecting random species from a distribution of species abundances using analytical
 * solutions from Vallade and Houchmandzadeh (2003) and Alonso and McKane (2004).
 *
 * Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
 */

#ifndef ANALYICAL_SPECIES_ABUNDANCES_H
#define ANALYICAL_SPECIES_ABUNDANCES_H

#include "SpeciesAbundancesHandler.h"
#include "neutral_analytical.h"
#include "RNGController.h"

namespace na = neutral_analytical;
using namespace std;

class AnalyticalSpeciesAbundancesHandler : public virtual SpeciesAbundancesHandler
{
protected:
    unsigned long seen_no_individuals;
    // Store all previous species ids in a map of cumulative numbers of individuals for searching for ids
    map<unsigned long, unsigned long> ind_to_species;
public:

    /**
     * @brief Default constructor
     */
    AnalyticalSpeciesAbundancesHandler();

    /**
     * @brief Default destructor
     */
    ~AnalyticalSpeciesAbundancesHandler() override = default;

    /**
     * @brief Creates the SpeciesAbundancesHandler object
     * @param random the random number generator
     * @param metacommunity_size the number of individuals in the metacommunity
     * @param speciation_rate the speciation rate of the metacommunity
     * @param local_community_size: the number of individuals in the local community
     */
    void setup(shared_ptr<RNGController> random, const unsigned long &metacommunity_size,
                   const long double &speciation_rate,
                   const unsigned long &local_community_size) override;

    /**
     * @brief Generates the species abundances using the analytical approximation.
     */
    void generateSpeciesAbundances();

    /**
     * @brief Gets a randomly generated species identity.
     * @return the species identity
     */
    unsigned long getRandomSpeciesID() override;

    /**
     * @brief Picks out a random individual from previously-seen individuals.
     * @param individual_id the individual id number to pick
     * @return species id of the individual
     */
    unsigned long pickPreviousIndividual(const unsigned long &individual_id);

    /**
     * @brief Picks out a new individual/species id with a random species abundance.
     * @return the species id of the new individual
     */
    void addNewSpecies();

    /**
     * @brief Gets a random species abundance by sampling from the logarithmic distribution.
     * @note this produces the abundance of any given species, not the abundance of any given individual
     * @return the randomly generated abundance of a species
     */
    unsigned long getRandomAbundanceOfSpecies();
};

#endif //ANALYICAL_SPECIES_ABUNDANCES_H
