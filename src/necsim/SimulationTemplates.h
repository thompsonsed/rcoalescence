// This file is part of necsim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.
//
/**
 * @author Sam Thompson
 * @date 02/01/2018
 * @file SimulationTemplates.h
 * @brief Contains template function for running any class of simulation (including protracted simulations, spatial and
 * non-spatial simulations.
 *
 * Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */

#ifndef SIMULATIONTEMPLATES_H
#define SIMULATIONTEMPLATES_H

#include <string>
#include <sstream>
#include "Logging.h"
#include "custom_exceptions.h"

/**
 * @brief Gets the [2] element from the vector (which should contain the config file from command-line arguments
 * @param com_args the vector of command-line arguments
 * @return the string at the [2] position containing the path to the config file
 */
const string &getConfigFileFromCmdArgs(const vector<string> &com_args)
{
    if(com_args.size() != 3)
    {
        stringstream ss;
        ss << "Incorrect number of command-line arguments supplied. Should be 3, got " << com_args.size() << endl;
        throw FatalException(ss.str());
    }
    else
    {
        return com_args[2];
    }
}

/**
 * @brief Template class for running simulations from all Tree types.
 * @tparam T the class (either Tree, or a child of Tree) of the simulation
 * @param config_file the config file to read simulation parameters from
 */
template<class T>
void runMain(const string &config_file)
{
    // Create our tree object that contains the simulation
    T tree;
    tree.importSimulationVariables(config_file);
    // Setup the sim
    tree.setup();
    // Detect speciation rates to apply
    bool isComplete = tree.runSimulation();
    if(isComplete)
    {
        tree.applyMultipleRates();
    }
    writeInfo("*************************************************\n");
}

#endif //SIMULATIONTEMPLATES_H
