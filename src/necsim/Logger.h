// This file is part of necsim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.
/**
 * @author Sam Thompson
 * @file Logger.h
 * @brief Routines for writing to cout. Intended to be overloaded for pythonic versions with the logging module.
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */
#ifndef LOGGING_IMPORT_H
#define LOGGING_IMPORT_H

#include <string>
#include <iostream>
#include <cstdio>
#include <stdexcept>
#include <sstream>
#include "LogFile.h"
#include "cpl_custom_handler.h"

using namespace std;

/**
 * @brief Controls writing to console and files for informing of statuses and debugging.
 */
class Logger
{
protected:
#ifdef DEBUG
    LogFile logfile;
#endif //DEBUG
public:
    Logger() = default;

    virtual ~Logger() = default;

    /**
     * @brief Logs an information message.
     * @param message the message to write out
     */
    virtual void writeInfo(string message);

    /**
     * @brief Logs a warning message.
     * @param message the message to write out
     */
    virtual void writeWarning(string message);

    /**
     * @brief Logs an error message.
     * @param message the message to write out
     */
    virtual void writeError(string message);

    /**
     * @brief Logs a critical message.
     * @param message the message to write out
     */
    virtual void writeCritical(string message);

#ifdef DEBUG

    /**
     * @brief Logs the message at the desired level.
     * @param level the level of logging severity
     * @param message the message to pass out as a string
     */
    virtual void writeLog(const int &level, string message);

    /**
     * @brief Logs the message at the desired level.
     * @param level the level of logging severity
     * @param message the message to pass out as a stringstream
     */
    virtual void writeLog(const int &level, stringstream &message);

#endif // DEBUG
};

#endif // LOGGING_IMPORT_H
