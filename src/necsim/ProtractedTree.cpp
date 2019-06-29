// This file is part of necsim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.
//
/**
 * @file ProtractedTree.cpp
 * @brief Contains the ProtractedTree class for running simulations and outputting the phylogenetic trees using protracted
 * speciation.
 *
 * Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */

#include "ProtractedTree.h"

bool ProtractedTree::calcSpeciation(const long double &random_number,
                                    const long double &speciation_rate,
                                    const unsigned long &no_generations)
{
    if(generation < speciation_generation_min)
    {
        return false;
    }
    if(generation >= speciation_generation_max)
    {
        return true;
    }
    return checkSpeciation(random_number, speciation_rate, no_generations);
}

void ProtractedTree::speciateLineage(const unsigned long &data_position)
{
    (*data)[data_position].setSpec(0.0);
    if(speciation_generation_min >= (*data)[data_position].getGenRate() + (*data)[data_position].getGeneration())
    {
        (*data)[data_position].setGenerationRate(static_cast<unsigned long>(floor(speciation_generation_min)) + 2);
    }
#ifdef DEBUG
    if(generation < speciation_generation_min)
    {
        (*data)[data_position].logLineageInformation(50);
        throw FatalException("Speciation attempted before minimum speciation generation. Please report this bug.");
    }
#endif // DEBUG
    (*data)[data_position].speciate();
}

bool ProtractedTree::getProtracted()
{
    return true;
}

void ProtractedTree::setProtractedVariables(double speciation_gen_min_in, double speciation_gen_max_in)
{
    speciation_generation_min = speciation_gen_min_in;
    speciation_generation_max = speciation_gen_max_in;
}

string ProtractedTree::getProtractedVariables()
{
    stringstream ss;
    ss << speciation_generation_min << "\n" << speciation_generation_max << "\n";
    return ss.str();
}

double ProtractedTree::getProtractedGenerationMin()
{
    return speciation_generation_min;
}

double ProtractedTree::getProtractedGenerationMax()
{
    return speciation_generation_max;
}

string ProtractedTree::protractedVarsToString()
{
    string tmp = "1 , " + to_string(speciation_generation_min) + ", " + to_string(speciation_generation_max);
    return tmp;
}