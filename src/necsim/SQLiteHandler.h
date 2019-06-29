// This file is part of necsim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @file SQLiteHandler.h
 * @brief Contains wrappers for the sqlite3 database objects for handling proper initiation, sharing between objects and
 * destruction.
 *
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */

#ifndef SQLITEHANDLER_H
#define SQLITEHANDLER_H

#include <sqlite3.h>
#include <utility>
#include <memory>

#ifdef WIN_INSTALL
#include <windows.h>
#define sleep Sleep
#else
#include <unistd.h>
#endif // WIN_INSTALL


#include "file_system.h"
#include "custom_exceptions.h"

/**
 * @brief Wrapper for sqlite3_stmt which stores the last command issued and the pointer to the statement object.
 */
struct SQLStatement
{
public:
    std::string last_command;
    sqlite3_stmt *stmt;

    /**
     * @brief Default constructor.
     */
    SQLStatement();

    /**
     * @brief Steps the transaction forwards one step.
     * @note Waits 10 seconds for locks to be released before raising an error if there are problems.
     * @return the result code of the transaction step
     */
    int step();

    /**
     * @brief Destroys the previous bindings and resets the statement.
     */
    void clearAndReset();

};

/**
 * Handler for the SQLite connection, including proper opening and closing of the database object.
 */
class SQLiteHandler
{
protected:
    sqlite3 *database;
    std::string file_name;
    shared_ptr<SQLStatement> stmt;
public:
    /**
     * @brief Default constructor
     */
    SQLiteHandler();

    /**
     * @brief Default destructor.
     */
    ~SQLiteHandler();

    /**
     * @brief Opens a database connection to the specified file name. If the file name is ":memory:", instead opens a
     * connection to an in-memory database object.
     * @param file_name the name of the file to open (or ":memory:" for in-memory databases)
     */
    void open(const std::string &file_name);

    /**
     * @brief Closes the sqlite3 connection to the database.
     */
    void close();

    /**
     * @brief Gets the error message from the sqlite3 database operations.
     * @param rc the result code of the previous operation, to print out
     * @return string containing the result code and the error message
     */
    std::string getErrorMsg(int rc);

    /**
     * @brief Gets the error message from the database
     * @return the error message
     */
    std::string getErrorMsg();

    /**
     * @brief Copies the data from the provided SQLiteHander object to this database.
     * @param sqlite_handler the database containing data to copy
     */
    void backupFrom(SQLiteHandler &sqlite_handler);

    /**
     * @brief Prepares the given commmand within the statement object.
     * @param command the command to execute
     * @param stmt the statement to prepare within
     * @return pointer to the prepared statement
     */
    shared_ptr<SQLStatement> prepare(const std::string &command);

    /**
     * @brief Creates a new statement for the database handler.
     */
    void createStatement();

    /**
     * @brief Use the supplied statement object for the database.
     * @param stmt
     */
    void useStatement(shared_ptr<SQLStatement> stmt);

    /**
     * @brief Steps the prepared statement forwards and reports any errors.
     * @note stmt should have been opened from the same SQLiteHandler object using prepare()
     * @param stmt the statement to step forwards
     */
    void step();

    /**
     * @brief Finalises the sqlite statement and reports any errors.
     * @note stmt should have been opened from the same SQLiteHandler object using prepare()
     * @param stmt the statement to finalise
     */
    void finalise();

    /**
     * @brief Executes a command from the database and reports any errors.
     * @param command the command to execute within the database
     */
    void execute(const string &command);

    /**
     * @brief Starts a transaction from this database object.
     */
    void beginTransaction();

    /**
     * @brief Ends the transaction from this database object.
     */
    void endTransaction();

    /**
     * @brief Checks if the database is open.
     * @return true, if the database is not a nullptr.
     */
    bool isOpen();

    /**
     * @brief Checks if the database has the specified table.
     * @param table_name the table name to check for existence
     * @return true if the table exists
     */
    bool hasTable(const std::string &table_name);
};

#endif //SQLITEHANDLER_H
