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

#include <sstream>
#include <stdexcept>
#include "neutral_analytical.h"
#include "RNGController.h"

namespace neutral_analytical
{
    long double siMetacommunitySpeciesWithAbundance(const unsigned long &n, const unsigned long &metacommunity_size,
                                                    const long double &speciation_rate)
    {
        long double fundamental_biodiversity_number = calcFundamentalBiodiversityNumber(metacommunity_size, speciation_rate);
        long double term1 = fundamental_biodiversity_number / n;
        long double term2 = lgamma(metacommunity_size + 1) + lgamma(metacommunity_size + fundamental_biodiversity_number - n);
        long double term3 = lgamma(metacommunity_size + 1 - n) + lgamma(metacommunity_size + fundamental_biodiversity_number);
        return term1 * exp(term2 - term3);
    }

    long double calcFundamentalBiodiversityNumber(const unsigned long &community_size,
                                                  const long double &speciation_rate)
    {
        if(speciation_rate == 1.0)
        {
            return 0.000000000000001;
        }

        return community_size * speciation_rate / (1.0 - speciation_rate);
    }

    long double calcSpeciationRate(const long double &fundamental_biodiversity_number,
                                   const unsigned long &metacommunity_size)
    {
        return 1.0/(1.0 + ((metacommunity_size - 1)/fundamental_biodiversity_number));
    }

    long double siSpeciesRichnessDeprecated(const unsigned long &community_size, const long double &speciation_rate)
    {
        long double total_species = 0.0;
        long double fundamental_biodiversity_number = calcFundamentalBiodiversityNumber(community_size, speciation_rate);
        for(unsigned long i = 0; i < community_size - 1; i++)
        {

            total_species += fundamental_biodiversity_number / (fundamental_biodiversity_number + i);
        }
        return total_species;
    }

    long double siSpeciesRichness(const unsigned long &community_size, const long double &speciation_rate)
    {
        return community_size * speciation_rate * log(1.0 / speciation_rate) / (1 - speciation_rate);
    }

    std::map<unsigned long, long double> siSpeciesAbundanceCumulativeDistribution(const unsigned long &community_size,
                                                                                  const long double &speciation_rate)
    {
        std::map<unsigned long, long double> sad;
        for(unsigned long i = 1; i < community_size; i++)
        {
            sad[i] = siMetacommunitySpeciesWithAbundance(i, community_size, speciation_rate);
        }
        return sad;
    }
}