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
namespace necsim
{
    class ProtractedSpatialTree : public SpatialTree, public ProtractedTree
    {
    public:
        ProtractedSpatialTree() : Tree(), SpatialTree(), ProtractedTree()
        {
        }

        ProtractedSpatialTree(ProtractedSpatialTree &&other) noexcept : ProtractedSpatialTree()
        {
            *this = std::move(other);
        }


        ~ProtractedSpatialTree() override= default;



        ProtractedSpatialTree(const ProtractedSpatialTree &other) : ProtractedSpatialTree()
        {
            *this = other;
        };

        ProtractedSpatialTree &operator=(ProtractedSpatialTree other) noexcept
        {
            other.swap(*this);
            return *this;
        }

        void swap(ProtractedSpatialTree &other) noexcept
        {
            if(this != &other)
            {
                SpatialTree::swap(other);
                std::swap(speciation_generation_min, other.speciation_generation_min);
                std::swap(speciation_generation_max, other.speciation_generation_max);

            }
        }
    };


}
#endif //SPECIATIONCOUNTER_PROTRACTEDSPATIALTREE_H
