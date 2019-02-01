// This file is part of necsim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details
/**
 * @author Samuel Thompson
 * @file Logging.cpp
 * @brief Routines for writing strings to the global logger object.
 * @copyright <a href="https://opensource.org/licenses/MIT">MIT Licence.</a>
 */

#include "Logging.h"
#include "Logger.h"

void writeInfo(string message)
{
    logger->writeInfo(message);
}

void writeWarning(string message)
{
    logger->writeWarning(message);
}

void writeError(string message)
{
    logger->writeError(message);
}

void writeCritical(string message)
{
    logger->writeCritical(message);
}

#ifdef DEBUG
void writeLog(const int &level, string message)
{
    logger->writeLog(level, message);
}

void writeLog(const int &level, stringstream &message)
{
    logger->writeLog(level, message.str());
}

#endif //DEBUG