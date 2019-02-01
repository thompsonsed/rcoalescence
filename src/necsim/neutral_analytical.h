// This file is part of necsim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @author Samuel Thompson
 * @file neutral_analytical.cpp
 *
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 * @brief Contains a set of functions associated with providing analytical solutions for neutral models in ecology.
 *
 * Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
 */

#ifndef NECSIM_NEUTRAL_ANALYTICAL_H
#define NECSIM_NEUTRAL_ANALYTICAL_H

#include <cmath>
#include <vector>
#include <map>

namespace neutral_analytical
{
    /**
     * @brief Calculate the expected number of species with abundance n within a non-spatially explicit
     * neutral model.
     *
     * Implementation of the SAD approximation from Vallade and Houchmandzadeh (2003).
     *
     * @param n the number of individuals expected within the species
     * @param community_size the total number of individuals in the community
     * @param speciation_rate the speciation rate
     * @return the number of species with abundance n
     */
    long double nseMetacommunitySpeciesWithAbundance(const unsigned long &n, const unsigned long &community_size,
                                                     const long double &speciation_rate);

    /**
     * @brief Gets the fundamental biodiversity number for the system.
     *
     * @param community_size the number of individuals in the community
     * @param speciation_rate the speciation rate
     * @return the fundamental biodiversity number
     */
    long double getFundamentalBiodiversityNumber(const unsigned long &community_size,
                                                 const long double &speciation_rate);

    /**
     * @brief Gets the expected species richness for a non-spatially explicit neutral model.
     *
     * Using the sampling formula from Vallade and ﻿Houchmandzadeh (2003).
     *
     * @param community_size the number of individuals in the community
     * @param speciation_rate the speciation rate
     * @deprecated this version uses an old inefficient method
     * @return the number of species expected to exist
     */
    long double nseSpeciesRichnessDeprecated(const unsigned long &community_size, const long double &speciation_rate);

    /**
     * @brief Gets the expected species richness for a non-spatially explicit neutral model.
     *
     * Using an analytical approximation identified here by Sam Thompson (but likely an old published result).
     *
     * @param community_size the number of individuals in the community
     * @param speciation_rate the speciation rate
     * @deprecated this version uses an old inefficient method
     * @return the number of species expected to exist
     */
    long double nseSpeciesRichness(const unsigned long &community_size, const long double &speciation_rate);

    /**
     * @brief Get the species abundances probability distribution  for a set of community parameters within a
     * non-spatial neutral model of ecology.
     *
     * The abundances will be cumulative probabilities.
     *
     * @param community_size the number of individuals in the community
     * @param speciation_rate the speciation rate
     * @return
     */
    std::map<unsigned long, long double> nseSpeciesAbundanceCumulativeDistribution(const unsigned long &community_size,
                                                                                   const long double &speciation_rate);

}
#endif //NECSIM_NEUTRAL_ANALYTICAL_H
