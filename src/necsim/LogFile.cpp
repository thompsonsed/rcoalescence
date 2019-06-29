// This file is part of necsim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.
/**
 * @author Sam Thompson
 * @file LogFile.cpp
 * @brief Contains a class for logging to a logfile, including reporting level and timestamps.
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */

#include <sstream>
#include <iomanip>
#include <boost/filesystem.hpp>
#include "LogFile.h"
#include "file_system.h"
#include "custom_exceptions.h"

string getTime()
{
#ifdef DEBUG
    auto t = time(nullptr);
    auto tm = *localtime(&t);
    stringstream ss;
    ss << put_time(&tm, LOGNAME_FORMAT);
    return ss.str();
#endif
    return string("");
}

string getDefaultLogFile()
{
    time_t now = time(0);
    static char name[30];
    strftime(name, sizeof(name), LOGNAME_FORMAT, localtime(&now));
    boost::filesystem::path full_path(boost::filesystem::current_path());
    string out = full_path.string() + "/log/" + string(name) + ".log";
    return out;
}

void getUniqueFileName(string &basic_string)
{
    boost::filesystem::path file_path(basic_string);
    const string file_name = basic_string.substr(0, basic_string.find('.'));
    const string file_extension = basic_string.substr(basic_string.find('.'));
    unsigned long iterator = 0;
    while(boost::filesystem::exists(file_path))
    {
        if(iterator > 10000000)
        {
            basic_string = file_name + file_extension;
            throw FatalException("Could not create unique file name after 10000000 tries.");
        }
        basic_string = file_name + "_" + to_string(iterator) + file_extension; // NOLINT
        file_path = boost::filesystem::path(basic_string);
        iterator++;
    }
}

LogFile::LogFile()
{
    init(getDefaultLogFile());
}

LogFile::LogFile(string file_name_in)
{
    init(file_name_in);
}

LogFile::~LogFile()
{

    output_stream << getTime() << " LOGGING ENDED" << endl;
    output_stream.close();
}

void LogFile::init(string file_name_in)
{
    createParent(file_name_in);
    file_name = file_name_in;
    getUniqueFileName(file_name);
    output_stream.open(file_name);
    if(!output_stream)
    {
        // Throw runtime error to avoid problems of attempting to write py object out before the logger has been set.
        throw runtime_error("Could not create log file at " + file_name_in + ".");
    }
    levels_map[0] = "noneset";
    levels_map[10] = "debug";
    levels_map[20] = "info";
    levels_map[30] = "warning";
    levels_map[40] = "error";
    levels_map[50] = "critical";
    output_stream << getTime() << " LOGGING STARTED" << endl;
}

void LogFile::write(const int &level, string message)
{
    if(levels_map.count(level) == 0)
    {
        throw FatalException("Logging level must be one of 0, 10, 20, 30, 40 or 50.");
    }
    output_stream << getTime() << " ";
    replace(message.begin(), message.end(), '\n', ' ');
    output_stream << levels_map[level] << ": " << message << endl;

}

void LogFile::write(const int &level, stringstream &message)
{
    write(level, message.str());
}



