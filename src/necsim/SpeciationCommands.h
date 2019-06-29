//This file is part of necsim project which is released under MIT license.
//See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @author Samuel Thompson
 * @date 19/07/2017
 * @file SpeciationCommands.h
 *
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 * @brief Contains the SpeciationCommands class for performing calculations of the coalescence tree structure and generating
 *  the SQL database objects.
 * 
 * For use on the SQL database outputs of necsim v3.1+. It requires command line parameters and generates a data object from them.
 * Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
 */


#include <cstdio>
#include <memory>
#include "Community.h"
#include "TreeNode.h"
#include "SpecSimParameters.h"

/**
 * @brief Routines for parsing command-line arguments to apply speciation rates to a necsim output
 * post-simulation.
 */
class SpeciationCommands
{
private:
    // Contains all speciation current_metacommunity_parameters
    shared_ptr<SpecSimParameters> sp;
    // Command-line arguments for parsing
    vector<string> comargs;
    // number of command-line arguments
    int argc;

public:

    /**
     * @brief Default constructor for SpeciationCommands class.
     */
    SpeciationCommands() : sp(make_shared<SpecSimParameters>())
    {

    }

    /**
     * @brief Run the command line arguments check.
     * Writes arguments to the SpecSimParameters object
     * @param argc the number of arguments.
     * @param comargs a vector filled with the command line arguments
     */
    void parseArgs();

    /**
     * @brief Runs the main program including parsing command line arguments and running the main analyses.
     * @param argc the number of command line arguments
     * @param argv the array of command line arguments
     * @return
     */
    int applyFromComargs(int argc_in, char **argv);

};