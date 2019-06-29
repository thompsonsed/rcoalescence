// This file is part of necsim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @author Samuel Thompson
 * @file SimulatedSpeciesAbundancesHandler.cpp
 *
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 * @brief Class for repeatedly selecting random species from a distribution of species abundances.
 *
 * Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
 */

#include "SimulatedSpeciesAbundancesHandler.h"

SimulatedSpeciesAbundancesHandler::SimulatedSpeciesAbundancesHandler()
        : species_abundances(), species_richness_per_abundance(),
          cumulative_abundance_map(
                  make_shared<map<unsigned long, unsigned long>>()), total_species_number(0), number_of_individuals(0){}

unsigned long SimulatedSpeciesAbundancesHandler::getRandomSpeciesID()
{
    unsigned long random_abundance = getRandomAbundanceOfIndividual();
    // TOD move this to debug
    if(species_abundances.count(random_abundance) == 0)
    {
        stringstream ss;
        ss << "No species abundances found for " << random_abundance << ". Please report this bug." << endl;
        throw FatalException(ss.str());
    }
    if(species_richness_per_abundance[random_abundance] == 0)
    {
        stringstream ss;
        ss << "Richness is 0 for abundance of " << random_abundance << endl;
        throw FatalException(ss.str());
    }
    if(species_richness_per_abundance[random_abundance] == 1)
    {
        return species_abundances[random_abundance][0];
    }
    unsigned long random_species_index = random->i0(species_richness_per_abundance[random_abundance] - 1);
    if(random_species_index >= species_abundances[random_abundance].size())
    {
        stringstream ss;
        ss << "Random species index larger than species abundance size: " << random_species_index << ">=";
        ss << species_abundances[random_abundance].size() << endl;
        throw FatalException(ss.str());
    }
    return species_abundances[random_abundance][random_species_index];
}

void SimulatedSpeciesAbundancesHandler::setAbundanceList(
        const shared_ptr<map<unsigned long, unsigned long>> &abundance_list_in)
{
    shared_ptr<vector<unsigned long>> abundance_list = make_shared<vector<unsigned long>>();
    abundance_list->reserve(abundance_list_in->size());
    metacommunity_size = 0;
    for(const auto &item: *abundance_list_in)
    {
        abundance_list->push_back(item.second);
        metacommunity_size += item.second;
    }
    generateAbundanceTable(abundance_list);
    generateCumulativeAbundances(abundance_list);
}

void SimulatedSpeciesAbundancesHandler::setAbundanceList(shared_ptr<vector<unsigned long>> abundance_list_in)
{
    max_species_id = 0;
    number_of_individuals = 0;
    generateAbundanceTable(abundance_list_in);
    generateCumulativeAbundances(abundance_list_in);
}

void SimulatedSpeciesAbundancesHandler::generateAbundanceTable(shared_ptr<vector<unsigned long>> abundance_list)
{
    writeInfo("Generating abundance table...");
    max_species_id = 0;
    for(const auto &item: (*abundance_list))
    {
        max_species_id++;
        species_abundances[item].push_back(max_species_id);
    }
    for(const auto &item: species_abundances)
    {
        species_richness_per_abundance[item.first] = item.second.size();
    }
    writeInfo("done.\n");
}

void SimulatedSpeciesAbundancesHandler::generateCumulativeAbundances(shared_ptr<vector<unsigned long>> abundance_list)
{
    writeInfo("Generating cumulative abundances...");
    if(abundance_list->empty())
    {
        throw FatalException("Abundance list is empty - please report this bug.");
    }
    // Holds a sorted abundance list (as doubles).
    vector<unsigned long> temp_cumulative_abundances;
    temp_cumulative_abundances.resize(abundance_list->size(), 0);
    partial_sum(abundance_list->begin(), abundance_list->end(), temp_cumulative_abundances.begin());
#ifdef DEBUG
    unsigned long total_sum = accumulate(abundance_list->begin(), abundance_list->end(), (unsigned long) 0);
    if(temp_cumulative_abundances.empty())
    {
        throw FatalException("Temporary cumulative abundance list is empty - please report this bug.");
    }
    if(total_sum == 0)
    {
        throw FatalException("Total sum of abundances is 0 - please report this bug.");
    }
#endif // DEBUG

    for(unsigned long i = 0; i < temp_cumulative_abundances.size(); i++)
    {
        (*cumulative_abundance_map)[temp_cumulative_abundances[i]] = (*abundance_list)[i];
    }
    if(cumulative_abundance_map->begin()->second == 0)
    {
#ifdef DEBUG
        if(cumulative_abundance_map->begin()->first != 0)
        {
            stringstream ss;
            ss << "Species abundance generated for " << cumulative_abundance_map->begin()->first;
            ss << " has zero abundance. Please report this bug." << endl;
            throw FatalException(ss.str());
        }
#endif // DEBUG
        // Remove the first 0 (if there is one)
        cumulative_abundance_map->erase(cumulative_abundance_map->begin());
    }
#ifdef DEBUG
    if(cumulative_abundance_map->rbegin()->first != metacommunity_size)
    {
        stringstream ss;
        ss << "Last cumulative abundance value (" << cumulative_abundance_map->rbegin()->first;
        ss << ") is not equal to community size (" << metacommunity_size << "). Please report this bug." << endl;
        throw FatalException(ss.str());
    }
    if(metacommunity_size == 1 && cumulative_abundance_map->size() != 1)
    {
        stringstream ss;
        ss << "Community size is 1, but cumulative abundance map has " << cumulative_abundance_map->size();
        ss << " elements. Please report this bug." << endl;
        throw FatalException(ss.str());
    }
#endif // DEBUG
    writeInfo("done.\n");
}

unsigned long SimulatedSpeciesAbundancesHandler::getRandomAbundanceOfIndividual()
{
#ifdef DEBUG
    if(cumulative_abundance_map->size() == 1)
    {
        if(cumulative_abundance_map->begin()->second == 0)
        {
            throw FatalException("Only one abundance found in abundance list, and it is 0. Please report this bug.");
        }
        return cumulative_abundance_map->begin()->second;
    }
    if(cumulative_abundance_map->empty())
    {
        throw FatalException("Cumulative abundance map is empty. Please report this bug.");
    }
    if(metacommunity_size == 1)
    {
        stringstream ss;
        ss << "Community size is 1 with cumulative abundance map: " << endl;
        for(const auto &item: *cumulative_abundance_map)
        {
            ss << item.first << ": " << item.second << endl;
        }
        ss << "Please report this bug." << endl;
        throw FatalException(ss.str());
    }
#endif // DEBUG
    return cumulative_abundance_map->upper_bound(random->i0(metacommunity_size - 1))->second;
}

