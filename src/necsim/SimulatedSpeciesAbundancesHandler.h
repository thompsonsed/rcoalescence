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
#include "NRrand.h"
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
    unsigned long getRandomAbundanceOfIndividual() override;

    /**
     * @brief Gets the species richness of a particular abundance class.
     * @param abundance the abundance class of the species
     * @return the number of species with that abundance
     */
    unsigned long getSpeciesRichnessOfAbundance(const unsigned long &abundance) override;
//
//	/**
//	 * @brief Generates an actual species abundance distribution from the probability distribution.
//	 *
//	 *
//	 * This process involves randomly adding species to their abundance class according to the probability distribution,
//	 * and is not guaranteed to generate the exact number of individuals as required by the metacommunity (but this
//	 * should be a good approximation for larger community sizes).
//	 */
//	void generateBiodiversity()
//	{
//		makeProbabilityDistributionCumulative();
//		species_abundances.clear();
//		while(number_of_individuals < community_size)
//		{
//			// Randomly add a species
//			double random_number = random_no->d01();
//			auto random_abundance_pointer = lower_bound(abundance_prob_dist->begin(), abundance_prob_dist->end(),
//														random_number);
//			unsigned long random_abundance = random_abundance_pointer - abundance_prob_dist->begin();
//			species_abundances.emplace_back(random_abundance);
//			number_of_individuals += random_abundance;
//		}
//
//	}
//
//	/**
//	 * @brief Counts the number of species expected to exist according to the probability distribution. Stores the
//	 * result in total_species_number.
//	 */
//	void countTotalSpecies()
//	{
//		total_species_number = 0.0;
//		if(abundance_prob_dist->empty())
//		{
//			throw FatalException("The species abundance distribution has no species. Please report this bug.");
//		}
//		(*abundance_prob_dist)[0] = 0.0;
//		for(const auto &item : *abundance_prob_dist)
//		{
//			total_species_number += item;
//		}
//	}
//
//	/**
//	 * @brief Makes the probability distribution cumulative, if it has not already been made so.
//	 */
//	void makeProbabilityDistributionCumulative()
//	{
//		countTotalSpecies();
//		if(!doubleCompare(total_species_number, 1.0, 0.00000000001))
//		{
//			for(unsigned long i = 1; i < abundance_prob_dist->size(); i++)
//			{
//				(*abundance_prob_dist)[i] += (*abundance_prob_dist)[i - 1];
//			}
//			for(unsigned long i = 1; i < abundance_prob_dist->size(); i++)
//			{
//				(*abundance_prob_dist)[i] /= total_species_number;
//			}
//			(*abundance_prob_dist)[0] = 0.0;
//		}
//	}
//
//	/**
//	 * @brief Makes the species abundances cumulative.
//	 *
//	 * Assumes that number_of_individuals holds the total number of individuals in the community.
//	 */
//	void makeSpeciesAbundancesCumulative()
//	{
//		cumulative_abundances->reserve(species_abundances.size());
//		partial_sum(species_abundances.begin(), species_abundances.end(), cumulative_abundances->begin(),
//					addAndRescale);
//	}
//
//	/**
//	 * @brief Gets a random species ID
//	 * @return a random species identity
//	 */
//	unsigned long getRandomSpeciesID()
//	{
//		// Select a random species
//		double random_number = random_no->d01();
//		auto species_id_pointer = lower_bound(species_abundances.begin(), species_abundances.end(), random_number);
//		return species_id_pointer - species_abundances.begin();
//	}

};

#endif //SIMULATED_SPECIES_ABUNDANCES_H
