//This file is part of necsim project which is released under MIT license.
//See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.
// 
/**
 * @author Sam Thompson
 * @file Setup.cpp
 * @brief Contains the command line parsing and setup options for necsim.
 * 
 * Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */

#include "setup.h"
#include "Logging.h"

// the old stdout
int saved_stdout;

#ifdef PROFILE

ofstream csv_output;

void logToCsv(string place, time_t start, time_t end)
{
    if(!csv_output.good())
    {
        try
        {
            csv_output.open("csvout.csv");
        }
        catch(exception &e)
        {
            throw Fatal_Main_Exception("Csv logger output not good: " + e.what());
        }
    }
    csv_output << place << "," << start << "," << end << endl;
}

void closeCsv()
{
    if(csv_output.good())
    {
        csv_output.close();
    }
}
#endif

void runAsDefault(vector<string> &comargs)
{
    writeInfo("Setting default variables on small grid.\n");
    comargs.push_back("-f");
    comargs.push_back("1");
    comargs.push_back("10");
    comargs.push_back("10");
    comargs.push_back("null");
    comargs.push_back("150");
    comargs.push_back("150");
    comargs.push_back("25");
    comargs.push_back("25");
    comargs.push_back("null");
    comargs.push_back("2000");
    comargs.push_back("2000");
    comargs.push_back("500");
    comargs.push_back("500");
    comargs.push_back("100");
    comargs.push_back("Default/");
    comargs.push_back("0.000009");
    comargs.push_back("2");
    comargs.push_back("1");
    comargs.push_back("1");
    comargs.push_back("4");
    comargs.push_back("1");
    comargs.push_back("0");
    comargs.push_back("100");
    comargs.push_back("null");
    comargs.push_back("null");
    comargs.push_back("0.5");
    comargs.push_back("20.0");
    comargs.push_back("2.0");
    comargs.push_back("null");
    comargs.push_back("null");
    comargs.push_back("0.000009");
}

void runLarge(vector<string> &comargs)
{
    writeInfo("Setting default variables on large grid.\n");
    comargs.push_back("-f");
    comargs.push_back("1");
    comargs.push_back("500");
    comargs.push_back("500");
    comargs.push_back("null");
    comargs.push_back("500");
    comargs.push_back("500");
    comargs.push_back("0");
    comargs.push_back("0");
    comargs.push_back("null");
    comargs.push_back("100");
    comargs.push_back("100");
    comargs.push_back("2500");
    comargs.push_back("2500");
    comargs.push_back("100");
    comargs.push_back("Default/");
    comargs.push_back("0.0001");
    comargs.push_back("2");
    comargs.push_back("10");
    comargs.push_back("1");
    comargs.push_back("3600");
    comargs.push_back("1");
    comargs.push_back("1");
    comargs.push_back("50000");
    comargs.push_back("null");
    comargs.push_back("null");
    comargs.push_back("0.5");
    comargs.push_back("20.0");
    comargs.push_back("2.0");
    comargs.push_back("null");
    comargs.push_back("null");
    comargs.push_back("0.001");
}

void runXL(vector<string> &comargs)
{
    writeInfo("Setting default variables on extra large grid.\n");
    comargs.push_back("-f");
    comargs.push_back("1");
    comargs.push_back("6000");
    comargs.push_back("6400");
    comargs.push_back("null");
    comargs.push_back("34000");
    comargs.push_back("28000");
    comargs.push_back("8800");
    comargs.push_back("14800");
    comargs.push_back("null");
    comargs.push_back("24000");
    comargs.push_back("20000");
    comargs.push_back("10320");
    comargs.push_back("8080");
    comargs.push_back("10");
    comargs.push_back("Default/");
    comargs.push_back("0.0000001");
    comargs.push_back("2");
    comargs.push_back("49");
    comargs.push_back("0.2");
    comargs.push_back("21600");
    comargs.push_back("1");
    comargs.push_back("3");
    comargs.push_back("600");
    comargs.push_back("null");
    comargs.push_back("null");
    comargs.push_back("0");
    comargs.push_back("2.2");
    comargs.push_back("1.0");
    comargs.push_back("null");
    comargs.push_back("null");
    comargs.push_back("0.000009");
}

void removeComOption(unsigned long &argc, vector<string> &comargs)
{
    // stupidly long species_id_list of possible arguments, but can't think of a better way to check this.
    if(comargs[1] == "-d" || comargs[1] == "-D" || comargs[1] == "-dl" || comargs[1] == "-dL" || comargs[1] == " -Dl" ||
       comargs[1] == "-DL" ||
       comargs[1] == "-dx" || comargs[1] == "-dX" || comargs[1] == "-DX" || comargs[1] == " -Dx" ||
       comargs[1] == "-c" || comargs[1] == "-C" ||
       comargs[1] == "-config" || comargs[1] == "-Config" || comargs[1] == "-f" || comargs[1] == "-h" ||
       comargs[1] == "-H" || comargs[1] == "-F)")
    {
        comargs.erase(comargs.begin() + 1);
        argc--;
    }
    return;
}




