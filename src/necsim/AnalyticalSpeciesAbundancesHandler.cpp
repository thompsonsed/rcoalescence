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

void AnalyticalSpeciesAbundancesHandler::setup(shared_ptr<NRrand> random, const unsigned long &community_size,
                                               const long double &speciation_rate)
{
    SpeciesAbundancesHandler::setup(random, community_size, speciation_rate);
    generateSpeciesAbundances();
}

void AnalyticalSpeciesAbundancesHandler::generateSpeciesAbundances()
{
    writeInfo("burning in species abundance...");
    while(seen_no_individuals < community_size)
    {
        addNewSpecies();
    }
    if(seen_no_individuals != community_size)
    {
        stringstream ss;
        ss << "Seen number of individuals (" << seen_no_individuals << ") does not match community size (";
        ss << community_size << ") - please report this bug" << endl;
        throw FatalException(ss.str());
    }
    writeInfo("done.\n");

}

unsigned long AnalyticalSpeciesAbundancesHandler::getRandomSpeciesID()
{
    // Either choose from previously seen individuals, or pick out a new individual of a new species.
    auto individual_id = random->i0(community_size - 1);
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

    unsigned long new_abundance = 0;
    do
    {
        new_abundance = getRandomAbundanceOfSpecies();
    }
    while(new_abundance > community_size - seen_no_individuals);
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
}

unsigned long AnalyticalSpeciesAbundancesHandler::getRandomAbundanceOfSpecies()
{
    // First generate a random abundance class
    return static_cast<unsigned long>(max(static_cast<double>(
                                                  min(random->randomLogarithmic(1.0 - speciation_rate),
                                                      community_size)), 1.0));
}


