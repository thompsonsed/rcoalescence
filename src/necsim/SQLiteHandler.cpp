// This file is part of necsim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @file SQLiteHandler.cpp
 * @brief Contains wrappers for the sqlite3 database objects for handling proper initiation, sharing between objects and
 * destruction.
 *
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */

#include "SQLiteHandler.h"

SQLStatement::SQLStatement() : last_command(), stmt(nullptr)
{

}

int SQLStatement::step()
{
    int rc = sqlite3_step(stmt);
    time_t start_check, end_check;
    time(&start_check);
    time(&end_check);
    // Note that a 10 second threshold is used to wait for output before raising an error.
    while(rc != SQLITE_DONE && (end_check - start_check) < 10 && rc != SQLITE_OK && rc != SQLITE_ROW &&
          (rc == SQLITE_BUSY || rc == SQLITE_LOCKED))
    {
        rc = sqlite3_step(stmt);
        sleep(1);
        time(&end_check);
    }
    return rc;
}

void SQLStatement::clearAndReset()
{
    sqlite3_clear_bindings(stmt);
    sqlite3_reset(stmt);
}

SQLiteHandler::SQLiteHandler() : database(nullptr), file_name(), stmt(nullptr)
{

}

SQLiteHandler::~SQLiteHandler()
{
    close();
}

void SQLiteHandler::open(const std::string &file_name)
{
    openSQLiteDatabase(file_name, database);
    this->file_name = file_name;
}

void SQLiteHandler::close()
{
    sqlite3_close_v2(database);
    database = nullptr;
}

std::string SQLiteHandler::getErrorMsg(int rc)
{
    std::stringstream ss;
    ss << "Error code: " << rc << ": " << getErrorMsg();
    return ss.str();
}

std::string SQLiteHandler::getErrorMsg()
{
    return sqlite3_errmsg(database);
}

void SQLiteHandler::backupFrom(SQLiteHandler &sqlite_handler)
{
    sqlite3_backup *backupdb = sqlite3_backup_init(database, "main", sqlite_handler.database, "main");
    int rc = sqlite3_backup_step(backupdb, -1);
    int counter = 0;
    while(counter < 10 && (rc == SQLITE_BUSY || rc == SQLITE_LOCKED))
    {
        counter++;
        rc = sqlite3_backup_step(backupdb, -1);
        sleep(1);
    }
    if(rc != SQLITE_OK && rc != SQLITE_DONE)
    {
        std::stringstream ss;
        ss << "Database backup cannot be completed to " << file_name << " from " << sqlite_handler.file_name;
        ss << ". Check write access on output folder." << endl << getErrorMsg(rc);
        throw FatalException(ss.str());
    }
    rc = sqlite3_backup_finish(backupdb);
    //			os << "rc: " << rc << endl;
    if(rc != SQLITE_DONE && rc != SQLITE_OK)
    {
        std::stringstream ss;
        ss << "Database backup cannot be finished to " << file_name << " from " << sqlite_handler.file_name;
        ss << ". Check write access on output folder." << endl << getErrorMsg(rc);
        throw FatalException(ss.str());
    }
}

shared_ptr<SQLStatement> SQLiteHandler::prepare(const std::string &command)
{
    createStatement();
    int rc = sqlite3_prepare_v2(database, command.c_str(), static_cast<int>(strlen(command.c_str())), &stmt->stmt,
                                nullptr);
    if(rc != SQLITE_OK && rc != SQLITE_DONE)
    {
        std::stringstream ss;
        ss << "Could not prepare statement: " << command << endl << getErrorMsg(rc);
        throw FatalException(ss.str());
    }
    stmt->last_command = command;
    return stmt;
}

void SQLiteHandler::createStatement()
{
    stmt = make_shared<SQLStatement>();
}

void SQLiteHandler::useStatement(shared_ptr<SQLStatement> stmt)
{
    this->stmt = stmt;
}

void SQLiteHandler::step()
{
    int rc = stmt->step();
    if(rc != SQLITE_OK && rc != SQLITE_DONE && rc != SQLITE_ROW)
    {
        std::stringstream ss;
        ss << "Could not step statement: " << stmt->last_command << endl;
        ss << getErrorMsg(rc);
        throw FatalException(ss.str());
    }
}

void SQLiteHandler::finalise()
{
    int rc = sqlite3_finalize(stmt->stmt);
    if(rc != SQLITE_OK && rc != SQLITE_DONE)
    {
        std::stringstream ss;
        ss << "Could not finalise statement: " << getErrorMsg(rc);
        throw FatalException(ss.str());
    }
}

void SQLiteHandler::execute(const string &command)
{
    int rc = sqlite3_exec(database, command.c_str(), nullptr, nullptr, nullptr);
    if(rc != SQLITE_OK && rc != SQLITE_DONE)
    {
        std::stringstream ss;
        ss << "Could not execute statement: " << command << endl;
        ss << getErrorMsg(rc);
        throw FatalException(ss.str());
    }
}

void SQLiteHandler::beginTransaction()
{
    execute("BEGIN TRANSACTION;");
}

void SQLiteHandler::endTransaction()
{
    execute("END TRANSACTION;");
}

bool SQLiteHandler::isOpen()
{
    return database != nullptr;
}

bool SQLiteHandler::hasTable(const std::string &table_name)
{
    string call1 = "select count(type) from sqlite_master where type='table' and name='" + table_name + "'";
    prepare(call1);
    step();
    auto has_table = static_cast<bool>(sqlite3_column_int(stmt->stmt, 0));
    finalise();
    return has_table;
}
