// This file is part of necsim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.
//
/**
 * @file ProtractedTree.h
 * @brief Contains the ProtractedTree class for running simulations and outputting the phylogenetic trees using
 * protracted speciation.
 *
 * Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 *
 */
#include <vector>
#include <string>

#include "SpatialTree.h"

#ifndef PROTRACTED_SPATIAL_TREE_H
#define PROTRACTED_SPATIAL_TREE_H

/**
 * @brief Contains the protracted tree class, for running simulations with procated speciation.
 */
class ProtractedTree : public virtual Tree
{
private:
    // Variables for the protracted speciation variables
    // The number of generations a lineage must exist before speciating.
    // Speciation is therefore not allowed before this time.
    // If this value is 0, it has no effect.
    double speciation_generation_min;
    // The number of generations a lineage can exist before speciating.
    // All remaining lineages are speciated at this time.
    // If this value is 0, it has no effect.
    double speciation_generation_max;
public:

    ProtractedTree() : Tree(), speciation_generation_min(0.0), speciation_generation_max(0.0)
    {
        bIsProtracted = true;
    }

    /**
     * @brief Calculates the speciation probability from the random number, speciation rate and number of generations a
     * lineage has existed for.
     * @param random_number the generated random number from 0-1
     * @param speciation_rate the speciation rate to be applied
     * @param no_generations the number of generations a lineage has existed for
     * @return if true, speciation has occured
     */
    bool calcSpeciation(const long double &random_number,
                        const long double &speciation_rate,
                        const unsigned long &no_generations) override;

    /**
     * @brief Performs the actual speciation.
     * Includes handling of speciated lineages under protracted conditions.
     * @param data_position the position in the array of TreeNodes for this lineage
     */
    void speciateLineage(const unsigned long &data_position) override;

    /**
     * @brief Gets the protractedness of the simulation.
     * Overridden by protracted child classes.
     * @return
     */
    bool getProtracted() override;

    /**
     * @brief Sets the protracted variables
     * @param speciation_gen_min the minimum number of generations to have passed before speciation is allowed
     * @param speciation_gen_max the maximum number of generations a lineage can exist for before it is speciated.
     */
    void setProtractedVariables(double speciation_gen_min, double speciation_gen_max) override;

    /**
     * @brief Gets the protracted variables and returns them as a single, newline separated string.

     * @return string containing the protracted variables, separated by newlines.
     */
    string getProtractedVariables() override;

    /**
     * @brief Gets the minimum number of generations a lineage must exist.
     *
     * @return double the number of generations a lineage must exist
     */
    double getProtractedGenerationMin() override;

    /**
     * @brief Gets the maximum number of generations a lineage can exist.
     *
     * @return double the number of generations a lineage must exist
     */
    double getProtractedGenerationMax() override;

    /**
     * @brief Outputs the protracted variables to a string.
     *
     * This function is intended to be overridden by derived classes.
     * It is intended the output is used for writing to SQL databases.
     *
     * @return string containing a list of the protracted speciation variables.
     */
    string protractedVarsToString() override;

    /**
     * @brief Applies the given speciation rate to the tree.
     *
     * @note Currently this just copies code from the version in tree, which is not ideal, but this avoids
     * creating an extra function.
     *
     * @param sr the required speciation rate
     */
    void applySpecRate(double sr, double t);
};

#endif // PROTRACTED_SPATIAL_TREE_H
