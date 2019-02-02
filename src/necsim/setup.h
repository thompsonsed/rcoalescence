//This file is part of necsim project which is released under MIT license.
//See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.
// 
/**
 * @author Sam Thompson
 * @file Setup.h
 * @brief Contains declarations for the command line parsing and setup options for necsim.
 * 
 * Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
 * @copyright <a href="https://opensource.org/licenses/MIT">MIT Licence.</a>
 */

#ifndef SETUP
#define SETUP

#include <string>
#include <vector>

#ifndef WIN_INSTALL

#include <unistd.h>

#endif

#include <sstream>
#include <ctime>
#include <boost/filesystem.hpp>
#include <cstdio>
#include <ctime>
#include <iostream>
#include <iomanip>

// Forward declaring the global variables
// store the log file name for access anywhere.
using namespace std;

extern string log_name;
// the old stdout
extern int saved_stdout;
#ifdef DEBUG
/**
 * @brief Opens the log file for redirecting stdout.
 * 
 * @param append if true, appends to the existing log file.
 */
void openLogFile(bool append);
#endif

/**
* @brief Sets up the command-line arguments for default parameters.
* 
* This is intended for testing purposes only.
* @param comargs a vector of command-line arguments for putting the parameters into.
*/
void runAsDefault(vector<string> &comargs);

/**
* @brief Sets up the command-line arguments for larger-scale default parameters.
* 
* This is intended for testing purposes only.
* @param comargs a vector of command-line arguments for putting the parameters into.
*/
void runLarge(vector<string> &comargs);

/**
* @brief Sets up the command-line arguments for default very large scale parameters.
* 
* This is intended for testing purposes only.
* @param comargs a vector of command-line arguments for putting the parameters into.
*/
void runXL(vector<string> &comargs);

/**
 * @brief Removes the command line options supplied, leaving just a clean vector with the correct data in. 

 */
void removeComOption(unsigned long &argc, vector<string> &comargs);

#endif // SETUP
