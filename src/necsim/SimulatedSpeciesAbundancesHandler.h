// This file is part of necsim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @author Samuel Thompson
 * @file SimulatedSpeciesAbundancesHandler.h
 *
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 * @brief Class for repeatedly selecting random species from a distribution of species abundances.
 *
 * Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
 */

#ifndef SIMULATED_SPECIES_ABUNDANCES_H
#define SIMULATED_SPECIES_ABUNDANCES_H

#include <string>
#include <memory>
#include <map>
#include <vector>
#include <cmath>
#include <numeric>

#include "neutral_analytical.h"
#include "RNGController.h"
#include "custom_exceptions.h"
#include "double_comparison.h"
#include "SpeciesAbundancesHandler.h"

namespace na = neutral_analytical;

using namespace std;

/**
 * @brief Generates species identities from a set of species abundances.
 */
class SimulatedSpeciesAbundancesHandler : public virtual SpeciesAbundancesHandler
{
protected:
    // Maps abundance values to a vector containing species ids
    map<unsigned long, vector<unsigned long>> species_abundances;
    // Maps abundance values to the maximum number of species expected to be contained.
    map<unsigned long, unsigned long> species_richness_per_abundance;
    // Maps cumulative probabilities of choosing each abundance to abundance values
    shared_ptr<map<unsigned long, unsigned long>> cumulative_abundance_map;
    // Total species number
    double total_species_number;
    unsigned long number_of_individuals;

public:
    /**
     * @brief Default constructor.
     */
    SimulatedSpeciesAbundancesHandler();

    /**
     * @brief Default destructor
     */
    ~SimulatedSpeciesAbundancesHandler() override = default;

    unsigned long getRandomSpeciesID() override;

    /**
     * @brief Sets the abundance list.
     * @param abundance_list_in map of species ids to abundances
     */
    void setAbundanceList(const shared_ptr<map<unsigned long, unsigned long>> &abundance_list_in) override;

    /**
     * @brief Sets the abundance list
     * @param abundance_list_in vector of species abundances
     */
    void setAbundanceList(shared_ptr<vector<unsigned long>> abundance_list_in) override;

    /**
     * @brief Generates the species abundance hash tables for efficiently storing the species identities.
     * @param abundance_list vector of species abundances
     */
    void generateAbundanceTable(shared_ptr<vector<unsigned long>> abundance_list);

    /**
     * @brief Generates the cumulative abundances, rescaled from 0-1, for randomly choosing an abundance using binary
     * search
     * @param abundance_list vector of species abundances
     */
    void generateCumulativeAbundances(shared_ptr<vector<unsigned long>> abundance_list);

    /**
     * @brief Gets a random species abundance.
     * @return the randomly generated abundance
     */
    unsigned long getRandomAbundanceOfIndividual();
};

#endif //SIMULATED_SPECIES_ABUNDANCES_H
