// This file is part of necsim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @author Samuel Thompson
 * @file SimulateDispersal.h
 * @brief Contains the ability to simulate a given dispersal kernel on a specified density map, outputting the effect 
 * dispersal distance distribution to an SQL file after n number of dispersal events (specified by the user).
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */

#ifndef DISPERSAL_TEST
#define DISPERSAL_TEST
#ifndef PYTHON_COMPILE
#define PYTHON_COMPILE
#endif

#include<string>
#include <cstdio>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include <sqlite3.h>
#include <set>
#include <mutex>
#include "Landscape.h"
#include "DispersalCoordinator.h"
#include "RNGController.h"
#include "Cell.h"
#include "DataMask.h"
#include "SQLiteHandler.h"

namespace necsim
{
    /**
     * @class SimulateDispersal
     * @brief Contains routines for importing a density map file, running a dispersal kernel n times on a landscape and record the
     * dispersal distances.
     */
    class SimulateDispersal
    {
    protected:
        // The density map object
        shared_ptr<Landscape> density_landscape;
        // The samplemask object
        DataMask data_mask;
        // Dispersal coordinator
        DispersalCoordinator dispersal_coordinator{};
        // Stores all key simulation current_metacommunity_parameters for the Landscape object
        shared_ptr<SimParameters> simParameters;
        // The random number generator object
        shared_ptr<RNGController> random;
        // The random number seed
        unsigned long seed;
        // The sqlite3 database object for storing outputs
        SQLiteHandler database;
        // Vector for storing pairs of dispersal distances to parameter references
        vector<std::tuple<unsigned long, Cell, double>> distances;
        // Maps distances to parameter references
        std::map<unsigned long, unsigned long> parameter_references;
        // Vector for storing the cells (for randomly choosing from)
        vector<Cell> cells;
        // The number of repeats to run the dispersal loop for
        unsigned long num_repeats;
        // The number of num_steps within each dispersal loop for the average distance travelled, which should be
        std::set<unsigned long> num_steps;
        // The number of threads launched to parallelise the distance simulation
        unsigned long num_workers;
        // generation counter
        double generation;
        // If true, sequentially selects dispersal probabilities, default is true
        bool is_sequential;
        // Reference number for this set of current_metacommunity_parameters in the database output
        unsigned long max_parameter_reference;

    public:

        SimulateDispersal() : density_landscape(make_shared<Landscape>()), data_mask(), dispersal_coordinator(),
                              simParameters(make_shared<SimParameters>()), random(make_shared<RNGController>()),
                              seed(0), database(), distances(), parameter_references(), cells(), num_repeats(0),
                              num_steps(), num_workers(), generation(0.0), is_sequential(false),
                              max_parameter_reference()
        {
        }

        ~SimulateDispersal() = default;

        SimulateDispersal(SimulateDispersal &&other) noexcept : SimulateDispersal()
        {
            *this = std::move(other);
        }

        SimulateDispersal(const SimulateDispersal &other) : SimulateDispersal()
        {
            density_landscape = other.density_landscape;
            data_mask = other.data_mask;
            dispersal_coordinator = other.dispersal_coordinator;
            simParameters = other.simParameters;
            random = other.random;
            seed = other.seed;
            database = other.database;
            distances = other.distances;
            parameter_references = other.parameter_references;
            cells = other.cells;
            num_repeats = other.num_repeats;
            num_steps = other.num_steps;
            num_workers = other.num_workers;
            generation = other.generation;
            is_sequential = other.is_sequential;
            max_parameter_reference = other.max_parameter_reference;
        };

        SimulateDispersal &operator=(SimulateDispersal other) noexcept
        {
            other.swap(*this);
            return *this;
        }

        void swap(SimulateDispersal &other) noexcept
        {
            if(this != &other)
            {
                std::swap(density_landscape, other.density_landscape);
                std::swap(data_mask, other.data_mask);
                other.dispersal_coordinator.swap(dispersal_coordinator);
                std::swap(simParameters, other.simParameters);
                std::swap(random, other.random);
                std::swap(seed, other.seed);
                std::swap(database, other.database);
                std::swap(distances, other.distances);
                std::swap(parameter_references, other.parameter_references);
                std::swap(cells, other.cells);
                std::swap(num_repeats, other.num_repeats);
                std::swap(num_steps, other.num_steps);
                std::swap(num_workers, other.num_workers);
                std::swap(generation, other.generation);
                std::swap(is_sequential, other.is_sequential);
                std::swap(max_parameter_reference, other.max_parameter_reference);
            }
        }

        /**
         * @brief Sets the is_sequential flag
         * @param bSequential if true, dispersal events are selected using the end point of the last dispersal
         * distance for the start of the next move event
         */
        void setSequential(bool bSequential);

        /**
         * @brief Sets the pointer to the simulation parameters object
         * @param sim_parameters pointer to the simulation parameters to use
         * @param p if true, writes the parameters out using writeInfo()
         */
        void setSimulationParameters(shared_ptr<SimParameters> sim_parameters, bool p = true);

        /**
         * @brief Import the maps from the simulation parameters.
         */
        void importMaps();

        /**
         * @brief Creates the map of steps to parameter references and initialises object sizes.
         */
        void setSizes();

        /**
         * @brief Sets the dispersal parameters in the DispersalCoordinator.
         */
        void setDispersalParameters();

        /**
         * @brief Sets the seed for the random number generator
         * @param s the seed
         */
        void setSeed(unsigned long s)
        {
            seed = s;
            random->wipeSeed();
            random->setSeed(s);
        }

        /**
         * @brief Sets the output database for writing results to
         * @param out_database path to the output database
         */
        void setOutputDatabase(string out_database);

        /**
         * @brief Sets the number of repeats to run the dispersal kernel for
         * @param n the number of repeats
         */
        void setNumberRepeats(unsigned long n);

        /**
         * @brief Sets the number of steps to run each repeat of the dispersal kernel for when recording mean distance
         * travelled
         * @param s a vector containing each step variable to record the distance at
         */
        void setNumberSteps(const vector<unsigned long> &s);

        /**
         * @brief Sets the number of threads launched to parallelise the distance simulation
         * @param n the number of workers
         */
        void setNumberWorkers(unsigned long n);

        /**
         * @brief Gets the maximum number of steps that is to be applied.
         * @return
         */
        unsigned long getMaxNumberSteps();

        /**
         * @brief Calculates the list of cells to choose randomly from
         */
        void storeCellList();

        /**
         * @brief Gets a random cell from the list of cells
         * @return a Cell object reference containing the x and y positions to choose from
         */
        const Cell &getRandomCell();

        /**
         * @brief Checks the density a given distance from the start point, calling the relevant landscape function.
         *
         * This also takes into account the rejection sampling of density based on the maximal density value from the map.
         *
         * @param this_cell Cell containing the x and y coordinates of the starting position
         */
        void getEndPoint(Cell &this_cell);

        /**
         * @brief Checks the density a given distance from the start point, calling the relevant landscape function.
         *
         * This also takes into account the rejection sampling of density based on the maximal density value from the map.
         *
         * @param this_cell Cell containing the x and y coordinates of the starting position
         * @param dispersal_coordinator Reference to the dispersal coordinator to use
         */
        void getEndPoint(Cell &this_cell, DispersalCoordinator &dispersal_coordinator);

        /**
         * @brief Runs the distance simulation edix-bidx times and reports the progress
         *
         * @tparam chooseRandomCells If true random walks will be chosen randomly from the cells vector, otherwise uses cells[bidx:edix]
         *
         * @param bidx First inclusive index into the cells vector of random walk origins to simulate in this worker
         * @param eidx Last exclusive index into the cells vector of random walk origins to simulate in this worker
         * @param num_repeats The number of repeats to average over for each cell
         * @param mutex The mutex to synchronise progress feedback to the user
         * @param finished The total number of cells simulated across all workers
         * @param dispersal_coordinator Reference to the dispersal corrdinator to use
         * @param generation Reference to the generation variable used byt the dispersal coordinator
         */
        template<bool chooseRandomCells = true> void runDistanceLoop(const unsigned long bidx,
                                                                     const unsigned long eidx,
                                                                     const unsigned long num_repeats,
                                                                     std::mutex &mutex,
                                                                     unsigned long &finished,
                                                                     DispersalCoordinator &dispersal_coordinator,
                                                                     double &generation);

        /**
         * @brief Runs the distance simulation eidx-bidx times on a separate worker and reports the progress
         *
         * @tparam chooseRandomCells If true random walks will be chosen randomly from the cells vector, otherwise uses cells[bidx:edix]
         *
         * @param seed Seed to use to initialse the per worker rng and dispersal coordinator
         * @param bidx First inclusive index into the cells vector of random walk origins to simulate in this worker
         * @param eidx Last exclusive index into the cells vector of random walk origins to simulate in this worker
         * @param num_repeats The number of repeats to average over for each cell
         * @param mutex The mutex to synchronise progress feedback to the user
         * @param finished The total number of cells simulated across all workers
         */
        template<bool chooseRandomCells = true> void runDistanceWorker(const unsigned long seed,
                                                                       const unsigned long bidx,
                                                                       const unsigned long eidx,
                                                                       const unsigned long num_repeats,
                                                                       std::mutex &mutex,
                                                                       unsigned long &finished);

        /**
         * @brief Simulates the dispersal kernel for the set parameters, storing the mean dispersal distance
         */
        void runMeanDispersalDistance();

        /**
         * @brief Simulates the dispersal kernel for the set parameters, storing the mean distance travelled.
         */
        void runMeanDistanceTravelled();

        /**
         * @brief Simulates the dispersal kernel for the set parameters on all habitable cells, storing the mean distance travelled.
         */
        void runAllDistanceTravelled();

        /**
         * @brief Simulates the dispersal kernel for the set parameters on the given sample cells, storing the mean distance travelled.
         *
         * @param samples Vector of cells to be sampled during the distance simulation
         */
        void runSampleDistanceTravelled(const vector<Cell> &samples);

        void runSampleDistanceTravelled(const vector<long> &sample_x, const vector<long> &sample_y);

        /**
         * @brief Writes the information about this repeat to the logger.
         */
        void writeRepeatInfo(unsigned long i);

        /**
         * @brief Writes out the distances to the SQL database.
         * @param table_name the name of the table to output to, either 'DISPERSAL_DISTANCE' or 'DISTANCES_TRAVELLED'
         */
        void writeDatabase(string table_name);

        /**
         * @brief Writes the simulation parameters to the output SQL database.
         * @param table_name the name of the table to output to, either 'DISPERSAL_DISTANCE' or 'DISTANCES_TRAVELLED'
         */
        void writeParameters(string table_name);

        /**
         * @brief Clears the parameters from the internal objects so that another dispersal simulation can be run, if
         * necessary.
         */
        void clearParameters();

        /**
         * @brief Gets the maximum parameter reference from the output SQL database and saves val + 1 to parameter_reference
         * Assumes that the database exists.
         *
         */
        void checkMaxParameterReference();

        /**
         * @brief Gets the maximum id number from the output SQL database and returns val + 1
         * Assumes that the database exists.
         * @note this function does not check for SQL injection attacks and should not be used with variable function names.
         * @param table_name: the name of the table to check for max(id) in
         * @returns the maximum id + 1 from the given table
         */
        unsigned long checkMaxIdNumber(string table_name);
    };
}
#endif
