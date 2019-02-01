// This file is part of necsim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @author Samuel Thompson
 * @file SpeciesAbundancesHandler.cpp
 *
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 * @brief Base class for storing and generating species abundances
 *
 * Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
 */

#include "SpeciesAbundancesHandler.h"
#include "custom_exceptions.h"

SpeciesAbundancesHandler::SpeciesAbundancesHandler() : random(make_shared<NRrand>()), max_species_id(0),
                                                       community_size(0), speciation_rate(0.0){}

void SpeciesAbundancesHandler::setup(shared_ptr<NRrand> random, const unsigned long &community_size,
                                     const long double &speciation_rate)
{
    SpeciesAbundancesHandler::random = std::move(random);
    SpeciesAbundancesHandler::community_size = community_size;
    SpeciesAbundancesHandler::speciation_rate = speciation_rate;
}

void SpeciesAbundancesHandler::setAbundanceList(const shared_ptr<map<unsigned long, unsigned long>> &abundance_list_in)
{

}

void SpeciesAbundancesHandler::setAbundanceList(shared_ptr<vector<unsigned long>> abundance_list_in)
{

}

unsigned long SpeciesAbundancesHandler::getRandomAbundanceOfIndividual()
{
    return 0;
}
