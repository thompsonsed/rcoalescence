// This file is part of necsim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @author Samuel Thompson
 * @date 19/07/2017
 * @file Filesystem.cpp
 *
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 * @brief Contains routines for checking files and folder exist, opening sqlite databases safely, with support for various
 * virtual filesystems, and checking parents of a file exist.
 *
 * Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
 */

#include <string>
#include <sstream>
#include <boost/filesystem.hpp>

#ifdef WIN_INSTALL
#include <windows.h>
#define sleep Sleep
#endif

#include "file_system.h"
#include "custom_exceptions.h"
#include "Logging.h"

void openSQLiteDatabase(const string &database_name, sqlite3 *&database)
{
    int rc;
    if(database_name == ":memory:")
    {
        rc = sqlite3_open(":memory:", &database);
        if(rc != SQLITE_OK && rc != SQLITE_DONE)
        {
            stringstream ss;
            ss << "Could not connect to in-memory database. Error: " << rc << endl;
            ss << " (" << sqlite3_errmsg(database) << ")" << endl;
            throw FatalException(ss.str());
        }
        return;
    }
#ifdef WIN_INSTALL
    rc = sqlite3_open_v2(database_name.c_str(), &database, SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE, "win32");
#else
    rc = sqlite3_open_v2(database_name.c_str(), &database, SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE, "unix-dotfile");
#endif
    if(rc != SQLITE_OK && rc != SQLITE_DONE)
    {
        int i = 0;
        while((rc != SQLITE_OK && rc != SQLITE_DONE) && i < 10)
        {
            i++;
            sleep(1);
            rc = sqlite3_open_v2(database_name.c_str(), &database, SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE,
                                 "unix-dotfile");
        }
        // Attempt different opening method if the first fails.
        int j = 0;
        while((rc != SQLITE_OK && rc != SQLITE_DONE) && j < 10)
        {
            j++;
            sleep(1);
            rc = sqlite3_open(database_name.c_str(), &database);
        }
        if(rc != SQLITE_OK && rc != SQLITE_DONE)
        {
            stringstream ss;
            ss << "ERROR_SQL_010: SQLite database file could not be opened. Check the folder exists and you "
                  "have write permissions. (REF1) Error code: "
               << rc << endl;
            ss << " Attempted call " << max(i, j) << " times" << endl;
            throw FatalException(ss.str());
        }
    }
}

void createParent(string file)
{
    // Boost < 1.59 support
    boost::filesystem::path file_path(file);
    auto parent = file_path.parent_path().string();
    if(parent.length() > 0)
    {
        std::string::iterator it = parent.end() - 1;
        if(*it == '/')
        {
            parent.erase(it);
        }
        boost::filesystem::path parent_path(parent);
        if(!boost::filesystem::exists(parent_path))
        {
            if(!parent_path.empty())
            {
                if(!boost::filesystem::create_directories(parent_path))
                {
                    throw runtime_error("Cannot create parent folder for " + file);
                }
            }
        }
    }
}

bool doesExist(string testfile)
{
    if(boost::filesystem::exists(testfile))
    {
        stringstream os;
        os << "\rChecking folder existance..." << testfile << " exists.               " << endl;
        writeInfo(os.str());
        return true;
    }
    else
    {
        throw runtime_error(string("ERROR_MAIN_008: FATAL. Input or output folder does not exist: " + testfile + "."));
    }
}

bool doesExistNull(string testfile)
{
    return testfile == "null" || testfile == "none" || doesExist(testfile);
}

unsigned long cantorPairing(const unsigned long &x1, const unsigned long &x2)
{
    return ((x1 + x2) * (x1 + x2 + 1) / 2) + x2;
}

unsigned long elegantPairing(const unsigned long &x1, const unsigned long &x2)
{
    if(x1 > x2)
    {
        return static_cast<unsigned long>(pow(x1, 2) + x1 + x2);
    }
    return static_cast<unsigned long>(pow(x2, 2) + x1);
}

vector<string> getCsvLineAndSplitIntoTokens(istream &str)
{
    vector<string> result;
    string line;
    getline(str, line);

    stringstream lineStream(line);
    string cell;

    while(getline(lineStream, cell, ','))
    {
        result.push_back(cell);
    }
    // This checks for a trailing comma with no data after it.
    if(!lineStream && cell.empty())
    {
        // If there was a trailing comma then add an empty element.
        result.emplace_back("");
    }
    return result;
}