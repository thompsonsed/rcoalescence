// This file is part of necsim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.
//
/**
 * @file ProtractedSpatialTree.h
 * @brief Contains the ProtractedSpatialTree class for running simulations and outputting the phylogenetic trees using
 * protracted speciation.
 *
 * Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 *
 */

#include "SpatialTree.h"
#include "ProtractedTree.h"

/**
 * @class ProtractedTree
 * @brief Contains the protracted tree class, for running simulations with procated speciation.
 */
#ifndef SPECIATIONCOUNTER_PROTRACTEDSPATIALTREE_H
#define SPECIATIONCOUNTER_PROTRACTEDSPATIALTREE_H

class ProtractedSpatialTree : public virtual SpatialTree, public virtual ProtractedTree
{
//    /**
//     * @brief Calculates the speciation probability from the random number, speciation rate and number of generations a
//     * lineage has existed for.
//     * @param random_number the generated random number from 0-1
//     * @param speciation_rate the speciation rate to be applied
//     * @param no_generations the number of generations a lineage has existed for
//     * @return if true, speciation has occured
//     */
//    bool calcSpeciation(const long double &random_number,
//                        const long double &speciation_rate,
//                        const unsigned long &no_generations) override
//    {
//        return ProtractedTree::calcSpeciation(random_number, speciation_rate, no_generations);
//    }

};

#endif //SPECIATIONCOUNTER_PROTRACTEDSPATIALTREE_H
