// This file is part of necsim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @author Samuel Thompson
 * @date 19/07/2017
 * @file SpeciationCommands.cpp
 *
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 * @brief Contains the ApplySpec class for performing calculations of the coalescence tree structure and generating
 *  the SQL database objects from the command-line.
 *
 * For use on the SQL database outputs of necsim v3.1+. It requires command line parameters and generates a data object
 * from them.
 * Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
 */

#include "SpeciationCommands.h"

void SpeciationCommands::parseArgs()
{
    bool bRunDefault = false;
    bool bInvalidArguments = false;
    bool bAskHelp = false;
    if(argc < 7)
    {
        if(argc == 1)
        {
            bInvalidArguments = true;
        }
        bInvalidArguments = true;
        if(argc == 2)
        {
            if((comargs[1]) == "-d")
            {
                bInvalidArguments = false;
                bRunDefault = true;
            }
            if(comargs[1] == "-h" || comargs[1] == "-help")
            {
                bInvalidArguments = false;
                bAskHelp = true;
            }
        }
        if(bInvalidArguments)
        {
            writeInfo("Incorrect number of arguments.");
            bInvalidArguments = true;
            if(argc == 1)
            {
                comargs.push_back("-e");
            }
            else
            {
                comargs[1] = "-e";
            }
        }
    }
    else
    {
        sp->samplemask = comargs[3];
        sp->filename = comargs[1];
        sp->times_file = comargs[4];
    }
    if(argc > 7)
    {
        sp->bMultiRun = true;
        int i = 6;
        while(i < argc)
        {
            sp->all_speciation_rates.insert(stof(comargs[i]));
            i++;
        }
    }
    else if(argc == 7 && !bInvalidArguments && !bAskHelp)
    {
        sp->bMultiRun = false;
        sp->all_speciation_rates.insert(stod(comargs[6]));
    }
    if(!bInvalidArguments && !bAskHelp && !bRunDefault)
    {
        if(comargs[2] == "true" || comargs[2] == "True" || comargs[2] == "T" || comargs[2] == "TRUE" ||
           comargs[2] == "t")
        {
            sp->use_spatial = true;
        }
        else
        {
            sp->use_spatial = false;
        }
        if(comargs[5] == "false" || comargs[5] == "False" || comargs[5] == "F" || comargs[5] == "FALSE" ||
           comargs[5] == "f")
        {
            sp->use_fragments = false;
        }
        else
        {
            if(comargs[5] == "true" || comargs[5] == "True" || comargs[5] == "T" || comargs[5] == "TRUE" ||
               comargs[5] == "t")
            {
                sp->fragment_config_file = "null";
            }
            else
            {
                sp->fragment_config_file = comargs[5];
            }
            sp->use_fragments = true;
        }
    }
    if(bInvalidArguments || bAskHelp)
    {
        stringstream os;
        os << "At least six command-line arguments are expected." << endl;
        os << "1 - Path to SQL database file." << endl;
        os << "2 - T/F of whether to record full spatial data." << endl;
        os << "3 - the sample mask to use (use null if no mask is to be used)" << endl;
        os << "4 - the file containing tempororal points of interest. If null, the present is used for all "
              "calculations."
           << endl;
        os << "5 - T/F of whether to calculate abundances for each rectangular fragment. Alternatively, provide a "
              "csv file with fragment data to be read."
           << endl;
        os << "6 - Speciation rate." << endl;
        os << "7 - onwards - Further speciation rates. [OPTIONAL]" << endl;
        os << "Would you like to run with the default paramenters?" << endl;
        os << "       (This requires a SQL database file at ../../Data/Coal_sim/Test_output/data_0_1.db)"
           << endl;
        os << "Enter Y/N: " << flush;
        writeInfo(os.str());
        string sDef;
        cin >> sDef;
        if(sDef == "Y" || sDef == "y")
        {
            bRunDefault = true;
        }
        else
        {
            bRunDefault = false;
            throw FatalException("Nothing to do!");
        }
    }
    if(comargs[1] == "-d" || bRunDefault)
    {
        sp->filename = "../../Data/Coal_sim/Test_output/data_0_1.db";
        sp->all_speciation_rates.insert(0.001);
        sp->samplemask = "null";
        sp->times_file = "null";
        sp->fragment_config_file = "null";
        sp->use_fragments = false;
        sp->use_spatial = true;
    }
}

int SpeciationCommands::applyFromComargs(int argc_in, char **argv)
{
    argc = argc_in;
    importArgs(static_cast<const unsigned int &>(argc), argv, comargs);
    parseArgs();
    Community tree_list;
    tree_list.apply(sp);
    return 0;
}

