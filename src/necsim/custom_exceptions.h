//This file is part of necsim project which is released under MIT license.
//See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.
// Author: Samuel Thompson
// Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
/**
 * @author Samuel Thompson
 * @file CustomExceptions.h
 * @brief Contains the various exceptions used by necsim.
 * 
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */

#ifndef CUSTOM_EXCEPTION_H
#define CUSTOM_EXCEPTION_H

#include <stdexcept>
#include <utility>
#include "Logging.h"

using namespace std;

/*!
 * @struct FatalException
 * @brief  This is called any time a fatal exception is called and the program is unwound and ended.
 */
struct FatalException : public runtime_error
{
    FatalException() : runtime_error("Fatal exception thrown at run time, quitting program. "){}

    /**
     * @brief Writes the message out to the logger if debug mode is enabled, otherwise just throws a runtime_error.
     * @param msg the message to write out and pass to runtime_error
     */
    explicit FatalException(string msg) : runtime_error(msg)
    {
#ifdef DEBUG
        writeLog(50, msg);
#endif //DEBUG
    }
};

/**
 * @brief A structure for all exceptions thrown within config processes.
 */
struct ConfigException : public FatalException
{
    ConfigException() : FatalException("Exception thrown at run time in config: "){};

    /**
     * @brief Generates a FatalException with the specified error message.
     * @param msg the message to pass to FatalException
     */
    explicit ConfigException(string msg) : FatalException(std::move(msg)){}
};

#endif // CUSTOM_EXCEPTION_H