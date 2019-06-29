// This file is part of necsim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @author Samuel Thompson
 * @file AnalyticalSpeciesAbundancesHandler.cpp
 *
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 * @brief Class for repeatedly selecting random species from a distribution of species abundances using analytical
 * solutions from Vallade and Houchmandzadeh (2003) and Alonso and McKane (2004).
 *
 * Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
 */

#include "AnalyticalSpeciesAbundancesHandler.h"
#include "custom_exceptions.h"

AnalyticalSpeciesAbundancesHandler::AnalyticalSpeciesAbundancesHandler() : seen_no_individuals(0), ind_to_species()
{

}

void AnalyticalSpeciesAbundancesHandler::setup(shared_ptr<RNGController> random,
                                               const unsigned long &metacommunity_size,
                                               const long double &speciation_rate,
                                               const unsigned long &local_community_size)
{
    SpeciesAbundancesHandler::setup(random, metacommunity_size, speciation_rate, local_community_size);
    generateSpeciesAbundances();
}

void AnalyticalSpeciesAbundancesHandler::generateSpeciesAbundances()
{
    writeInfo("Burning in species abundances...\n");
    auto expected_richness = std::max(
            static_cast<unsigned long>(neutral_analytical::siSpeciesRichness(metacommunity_size,
                                                                             speciation_rate)), (unsigned long) 1);
    // We use an approximation if the metacommunity richness is much larger than the local community size.
    if(expected_richness > 100 * local_community_size)
    {
        // Recalculate the fundamental biodiversity number to produce an identical result
        auto original_fbn = neutral_analytical::calcFundamentalBiodiversityNumber(metacommunity_size, speciation_rate);
        // Save the original metacommunity size and speciation rate.
        long double original_speciation_rate = speciation_rate;
        unsigned long original_metacommunity_size = metacommunity_size;
        // Reformat the metacommunity speciation rate to generate an equivalent metacommunity.
        metacommunity_size = 100 * local_community_size;
        speciation_rate = neutral_analytical::calcSpeciationRate(original_fbn, metacommunity_size);
        expected_richness = static_cast<unsigned long>(neutral_analytical::siSpeciesRichness(metacommunity_size,
                                                                                             speciation_rate));
        stringstream ss;
        ss << "\tRescaling large metacommunity size using fundamental biodiversity number..." << endl;
        ss << "\tNew metacommunity size: " << metacommunity_size << endl;
        ss << "\tNew speciation rate: " << speciation_rate << endl;
        ss << "\tGenerating abundances for " << expected_richness << " species" << endl;
        writeInfo(ss.str());
        for(unsigned long i = 0; i < expected_richness; i++)
        {
            addNewSpecies();
        }
        // Now replace the main variables
        metacommunity_size = original_metacommunity_size;
        speciation_rate = original_speciation_rate;
    }
    else
    {
        stringstream ss;
        ss << "\tGenerating abundances for " << expected_richness << " species" << endl;
        writeInfo(ss.str());
        // Otherwise just generate the full SAD - note that this only approximates the desired metacommunity size.
        for(unsigned long i = 0; i < expected_richness; i++)
        {
            addNewSpecies();
        }
    }
    // Make sure that we've seen at least as many individuals as in the local community.
    if(seen_no_individuals < local_community_size && metacommunity_size > local_community_size)
    {
        stringstream ss;
        ss << "Seen number of individuals (" << seen_no_individuals << ") is not more than local community size (";
        ss << local_community_size << ") - please report this bug" << endl;
        ss << "Metacommunity: " << endl;
        ss << "\tsize: " << metacommunity_size << endl;
        ss << "\tspeciation rate: " << speciation_rate << endl;
        ss << "Local community size: " << local_community_size << endl;
        throw FatalException(ss.str());
    }

}

unsigned long AnalyticalSpeciesAbundancesHandler::getRandomSpeciesID()
{
    // Select a random individual from the seen number of individuals
    auto individual_id = random->i0(seen_no_individuals - 1);
#ifdef DEBUG
    if(individual_id > seen_no_individuals)
    {
        stringstream ss;
        ss << "Random individual ID (" << individual_id << ") is greater than the number of seen indiviudals (";
        ss << seen_no_individuals << ")" << endl;
        throw FatalException(ss.str());
    }
    if(ind_to_species.empty())
    {
        throw FatalException(
                "No individuals have been seen yet, but an individual ID was generated. Please report this bug.");
    }
#endif // DEBUG
    return pickPreviousIndividual(individual_id);
}

unsigned long AnalyticalSpeciesAbundancesHandler::pickPreviousIndividual(const unsigned long &individual_id)
{
    return ind_to_species.upper_bound(individual_id)->second;
}

void AnalyticalSpeciesAbundancesHandler::addNewSpecies()
{
    max_species_id++;
    unsigned long new_abundance = getRandomAbundanceOfSpecies();
    unsigned long cumulative_abundance;
    if(seen_no_individuals == 0)
    {
        cumulative_abundance = new_abundance;
    }
    else
    {
        cumulative_abundance = new_abundance + ind_to_species.rbegin()->first;
    }
    ind_to_species[cumulative_abundance] = max_species_id;
    seen_no_individuals += new_abundance;
#ifdef DEBUG
    if(ind_to_species.rbegin()->first != seen_no_individuals)
    {
        stringstream ss;
        ss << "ind_to_species end does not equal seen no inds: " << ind_to_species.rbegin()->first << "!=";
        ss << seen_no_individuals << endl;
        throw FatalException(ss.str());
    }
    if(ind_to_species.rbegin()->second != max_species_id)
    {
        stringstream ss;
        ss << "Last species id has not been set correctly: " << ind_to_species.rbegin()->second << "!=";
        ss << max_species_id << endl;
        throw FatalException(ss.str());
    }
#endif //DEBUG
}

unsigned long AnalyticalSpeciesAbundancesHandler::getRandomAbundanceOfSpecies()
{
    // First generate a random abundance class
    return static_cast<unsigned long>(max(static_cast<double>(
                                                  min(random->randomLogarithmic(1.0 - speciation_rate),
                                                      metacommunity_size)), 1.0));
}



