// This file is part of necsim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details
/**
 * @author Samuel Thompson
 * @file Logging.h
 * @brief Routines for writing strings to the global logger object.
 * @copyright <a href="https://opensource.org/licenses/MIT">MIT Licence.</a>
 */

#ifndef NECSIM_LOGGING_H
#define NECSIM_LOGGING_H

#include <string>
#include "Logger.h"
// Global declaration of logger
/**
 * @brief Global object for logging
 */
extern Logger *logger;

/**
 * @brief Checks if the global logger object has been set up.
 */
bool loggerIsSetup();

/**
	 * @brief Writes to cout, or to info in logging module if being compiled with python
	 * @param message the message to write out
	 */
void writeInfo(string message);

/**
 * @brief Writes to cerr, or to warning in logging module if being compiled with python
 * @param message the message to write out
 */
void writeWarning(string message);

/**
 * @brief Writes to cerr, or to error in logging module if being compiled with python
 * @param message the message to write out
 */
void writeError(string message);

/**
 * @brief Writes to cerr, or to critical in logging module if being compiled with python
 * @param message the message to write out
 */
void writeCritical(string message);

#ifdef DEBUG

/**
 * @brief Calls the static logger object for logging out
 * @param level the level of logging severity
 * @param message the message to pass out as a string
 */
void writeLog(const int &level, string message);

/**
 * @brief Calls the static logger object for logging out
 * @param level the level of logging severity
 * @param message the message to pass out as a stringstream
 */
void writeLog(const int &level, stringstream &message);

#endif // DEBUG

#endif //MEANDISTANCEMODULE_LOGGING_H
