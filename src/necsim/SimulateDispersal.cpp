// This file is part of necsim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @author Samuel Thompson
 * @file SimulateDispersal.cpp
 * @brief Contains the ability to simulate a given dispersal kernel on a specified density map, outputting the effect
 * dispersal distance distribution to an SQL file after n number of dispersal events (specified by the user).
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */


#include <utility>

#include "SimulateDispersal.h"
#include "Logging.h"
#include "custom_exceptions.h"
#include "file_system.h"


void SimulateDispersal::setSequential(bool bSequential)
{
    is_sequential = bSequential;
}

void SimulateDispersal::setSimulationParameters(shared_ptr<SimParameters> sim_parameters, bool print)
{
    if(print)
    {
        writeInfo("********************************\nSetting simulation current_metacommunity_parameters...\n");
    }
    simParameters = std::move(sim_parameters);
    if(simParameters == nullptr)
    {
        throw FatalException("Simulation parameters are nullptr. Please report this bug.");
    }
    if(print)
    {
        simParameters->printSpatialVars();
    }
    simParameters->deme = 1.0;
    simParameters->deme_sample = 1.0;
}

void SimulateDispersal::importMaps()
{
    writeInfo("Starting map import...\n");
    if(simParameters == nullptr)
    {
        throw FatalException("Simulation current_metacommunity_parameters have not been set.");
    }
    density_landscape->setDims(simParameters);
    dispersal_coordinator.setMaps(density_landscape);
    setDispersalParameters();
    density_landscape->calcFineMap();
    density_landscape->calcCoarseMap();
    density_landscape->calcOffset();
    density_landscape->calcHistoricalFineMap();
    density_landscape->calcHistoricalCoarseMap();
    density_landscape->setLandscape(simParameters->landscape_type);
    density_landscape->recalculateHabitatMax();
    data_mask.importSampleMask(simParameters);
}

void SimulateDispersal::setSizes()
{
    checkMaxParameterReference();
    unsigned long i = max_parameter_reference;
    if(num_steps.empty())
    {
        num_steps.insert(1);
    }
    for(const auto &item : num_steps)
    {
        parameter_references.insert(make_pair(item, i));
        i++;
    }
    distances.clear();
    distances.resize(num_repeats * num_steps.size());

}

void SimulateDispersal::setDispersalParameters()
{
    dispersal_coordinator.setRandomNumber(random);
    dispersal_coordinator.setGenerationPtr(&generation);
    dispersal_coordinator.setDispersal(simParameters);

}

void SimulateDispersal::setOutputDatabase(string out_database)
{
    // Check the file is a database
    if(out_database.substr(out_database.length() - 3) != ".db")
    {
        throw FatalException("Output database is not a .db file, check file name.");
    }
    // Open our SQL connection to the database
    database.open(out_database);
    checkMaxParameterReference();
}

void SimulateDispersal::setNumberRepeats(unsigned long n)
{
    num_repeats = n;
}

void SimulateDispersal::setNumberSteps(const vector<unsigned long> &s)
{
    for(const auto item : s)
    {
        num_steps.insert(item);
    }
}

unsigned long SimulateDispersal::getMaxNumberSteps()
{
    unsigned long max_number_steps = 0;
    if(!num_steps.empty())
    {
        max_number_steps = *num_steps.rbegin();
    }
    else
    {
        throw FatalException("No steps have been set.");
    }
    return max_number_steps;
}

void SimulateDispersal::storeCellList()
{
    unsigned long total = 0;
    unsigned long cell_total = 0;
    // First count the number of density cells and pick a cell size
    for(unsigned long i = 0; i < simParameters->sample_y_size; i++)
    {
        for(unsigned long j = 0; j < simParameters->sample_x_size; j++)
        {
            if(data_mask.getVal(j, i, 0, 0))
            {
                cell_total++;
                total += density_landscape->getVal(j, i, 0, 0, 0.0);
            }
        }
    }
    writeInfo("Choosing from " + to_string(cell_total) + " cells of " + to_string(total) + " individuals.\n");
    cells.resize(total);
    unsigned long ref = 0;
    for(unsigned long i = 0; i < simParameters->sample_y_size; i++)
    {
        for(unsigned long j = 0; j < simParameters->sample_x_size; j++)
        {
            for(unsigned long k = 0; k < density_landscape->getVal(j, i, 0, 0, 0.0); k++)
            {
                cells[ref].x = j;
                cells[ref].y = i;
                ref++;
            }
        }
    }
}

const Cell &SimulateDispersal::getRandomCell()
{
    if(cells.size() == 0)
    {
        throw FatalException("No cells in landscape to simulate.");
    }
    auto index = static_cast<unsigned long>(floor(random->d01() * cells.size()));
    return cells[index];
}

void SimulateDispersal::getEndPoint(Cell &this_cell)
{
    Step tmp_step(this_cell);
    dispersal_coordinator.disperse(tmp_step);
    this_cell.x = tmp_step.oldx + tmp_step.oldxwrap * simParameters->sample_x_size;
    this_cell.y = tmp_step.oldy + tmp_step.oldywrap * simParameters->sample_y_size;
//	return (this->*getValFptr)(dist, angle, this_cell, end_cell);
}

void SimulateDispersal::runMeanDispersalDistance()
{
    writeInfo("Simulating dispersal " + to_string(num_repeats) + " times.\n");
    storeCellList();
    Cell this_cell{};
    this_cell = getRandomCell();
    // Set up the parameter reference
    setSizes();
    for(unsigned long i = 0; i < num_repeats; i++)
    {
        Cell start_cell{};
        if(!is_sequential)
        {
            // This takes into account rejection sampling based on density due to
            // setup process for the cell list
            this_cell = getRandomCell();
        }
        start_cell = this_cell;
        // Check the end point
        getEndPoint(this_cell);
        // Now store the output location
        distances[i] = make_pair(1, distanceBetweenCells(this_cell, start_cell));
    }
    writeInfo("Dispersal simulation complete.\n");
}

void SimulateDispersal::runMeanDistanceTravelled()
{
    stringstream ss;
    ss << "Simulating dispersal " << num_repeats << " times for (";
    // The maximum number of steps
    setSizes();
    unsigned long max_number_steps = getMaxNumberSteps();
    for(const auto &item : num_steps)
    {
        ss << item;
        if(item != max_number_steps)
        {
            ss << ", ";
        }
    }
    ss << ") generations.\n";
    writeInfo(ss.str());
    storeCellList();
    Cell this_cell{}, start_cell{};
    // Reference for the distances vector
    unsigned long dist_i = 0;
    for(unsigned long i = 0; i < num_repeats; i++)
    {
        writeRepeatInfo(i);
        // iterator for elements in the set.
        auto step_iterator = num_steps.begin();
        this_cell = getRandomCell();
        start_cell = this_cell;
        generation = 0.0;
        // Keep looping until we get a valid end point
        for(unsigned long j = 1; j <= max_number_steps; j++)
        {
            getEndPoint(this_cell);
            generation += 0.5;
            if(j == *step_iterator)
            {
                distances[dist_i] = make_pair(j, distanceBetweenCells(start_cell, this_cell));
                step_iterator++;
                dist_i++;
            }
        }
    }
    writeRepeatInfo(num_repeats);
    writeInfo("\nDispersal simulation complete.\n");
}

void SimulateDispersal::writeRepeatInfo(unsigned long i)
{
    stringstream os;
    os << "\rSimulating dispersal " << i << "/" << num_repeats;
    writeInfo(os.str());
}

void SimulateDispersal::writeDatabase(string table_name)
{
    if(database.isOpen())
    {
        if(table_name != "DISTANCES_TRAVELLED" && table_name != "DISPERSAL_DISTANCES")
        {
            string message = "Table name " + table_name;
            message += "  is not one of 'DISTANCES_TRAVELLED' or 'DISPERSAL_DISTANCES'.";
            throw FatalException(message);
        }
        // Write out the current_metacommunity_parameters
        checkMaxParameterReference();
        writeParameters(table_name);
        // Do the sql output
        // First create the table
        string create_table = "CREATE TABLE IF NOT EXISTS " + table_name + " (id INT PRIMARY KEY not null, ";
        create_table += " distance DOUBLE not null, parameter_reference INT NOT NULL);";
        database.execute(create_table);
        string insert_table = "INSERT INTO " + table_name + " (id, distance, parameter_reference) VALUES (?, ?, ?);";
        auto stmt = database.prepare(insert_table);
        database.beginTransaction();
        unsigned long max_id = checkMaxIdNumber(table_name);
        database.useStatement(stmt); // this could be cleaned up if checkMaxIDNumber comes before the insert statement.
        for(unsigned long i = 0; i < distances.size(); i++)
        {
            sqlite3_bind_int(stmt->stmt, 1, static_cast<int>(max_id + i));
            auto iter = parameter_references.find(distances[i].first);

#ifdef DEBUG
            if(iter == parameter_references.end())
            {
                throw FatalException("Cannot find parameter reference. Please report this bug.");
            }
#endif // DEBUG
            unsigned long reference = iter->second;
            if(reference > max_parameter_reference)
            {
                max_parameter_reference = reference;
            }
            sqlite3_bind_double(stmt->stmt, 2, distances[i].second);
            sqlite3_bind_int(stmt->stmt, 3, static_cast<int>(reference));
            int step = stmt->step();
            if(step != SQLITE_DONE)
            {
                stringstream ss;
                ss << "Could not insert into database." << endl;
                ss << database.getErrorMsg(step);
                throw FatalException(ss.str());
            }
            stmt->clearAndReset();
        }
        database.endTransaction();
        database.finalise();
    }
    else
    {
        throw FatalException("Database connection has not been opened, check programming.");
    }
    clearParameters();
}

void SimulateDispersal::writeParameters(string table_name)
{
    // Now add the current_metacommunity_parameters
    string create_table = "CREATE TABLE IF NOT EXISTS PARAMETERS (ref INT PRIMARY KEY not null,";
    create_table += "simulation_type TEXT not null, ";
    create_table += " sigma DOUBLE not null, tau DOUBLE not null, m_prob DOUBLE not null, cutoff DOUBLE NOT NULL,";
    create_table += "dispersal_method TEXT not null, map_file TEXT not null, seed INT NOT NULL, number_steps ";
    create_table += "INT NOT NULL, number_repeats INT NOT NULL);";
    database.execute(create_table);
    for(const auto &item : parameter_references)
    {
        string insert_table =
                "INSERT INTO PARAMETERS VALUES(" + to_string(item.second) + ", '" + table_name + "',";
        insert_table += to_string((long double) simParameters->sigma) + ",";
        insert_table +=
                to_string((long double) simParameters->tau) + ", " + to_string((long double) simParameters->m_prob);
        insert_table +=
                ", " + to_string((long double) simParameters->cutoff) + ", '" + simParameters->dispersal_method + "','";
        insert_table += simParameters->fine_map_file + "', " + to_string(seed) + ", " + to_string(item.first) + ", ";
        insert_table += to_string(num_repeats) + ");";
        database.execute(insert_table);
    }
}

void SimulateDispersal::clearParameters()
{
    distances.clear();
    parameter_references.clear();
    num_steps.clear();
}

void SimulateDispersal::checkMaxParameterReference()
{
    if(database.hasTable("PARAMETERS"))
    {
        string to_exec = "SELECT CASE WHEN COUNT(1) > 0 THEN MAX(ref) ELSE 0 END AS [Value] FROM PARAMETERS;";
        auto stmt = database.prepare(to_exec);
        database.step();
        max_parameter_reference = static_cast<unsigned long>(sqlite3_column_int(stmt->stmt, 0) + 1);
        // close the old statement
        database.finalise();
    }
    else
    {
        max_parameter_reference = 1;
    }

}

unsigned long SimulateDispersal::checkMaxIdNumber(string table_name)
{
    string to_exec = "SELECT CASE WHEN COUNT(1) > 0 THEN MAX(id) ELSE 0 END AS [Value] FROM " + table_name + ";";
    auto stmt = database.prepare(to_exec);
    database.step();
    auto max_id = static_cast<unsigned long>(sqlite3_column_int(stmt->stmt, 0) + 1);
    database.finalise();
    return max_id;
}


