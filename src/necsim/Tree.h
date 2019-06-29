// This file is part of necsim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @file Tree.h
 * @brief  Contains the main simulation object for spatially implicit coalescence simulations.
 *
 * Provides the basis for spatially explicit versions in SpatialTree, and protracted speciation versions in
 * ProtractedTree and ProtractedSpatialTree.
 *
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */
#ifndef TREE_H
#define TREE_H
//#ifndef sql_ram
//#define sql_ram
//#endif

#include <sqlite3.h>
#include <string>

#ifdef CXX14_SUPPORT
#include "memory.h"
#else

#include <memory>

#endif

#include "TreeNode.h"
#include "Matrix.h"
#include "SimParameters.h"
#include "RNGController.h"
#include "DataPoint.h"
#include "Community.h"
#include "file_system.h"
#include "custom_exceptions.h"
#include "Step.h"
#include "SQLiteHandler.h"

/**
 * @brief Main simulation class for performing a non-spatial neutral simulation and generating the phylogenetic tree of
 * the individuals.
 */
class Tree
{
protected:
    // storing the coalescence tree itself
    shared_ptr<vector<TreeNode>> data;
    // a reference for the last written point in data.
    unsigned long enddata;
    // Stores the command line current_metacommunity_parameters and parses the required information.
    shared_ptr<SimParameters> sim_parameters;
    // random number generator
    shared_ptr<RNGController> NR;
    // Storing the speciation rates for later reference.
    vector<long double> speciation_rates;
    // flag for having set the simulation seed.
    bool seeded;
    // random seed
    long long seed;
    // for file naming - good to know which job_type in a series is being executed here
    long long job_type;
    // The map file containing the times that we want to expand the model and record all lineages again.
    // If this is null, uses_temporal_sampling will be false and the vector will be empty.
    string times_file;
    vector<double> reference_times;
    // Set to true if we are recording at times other than the present day.
    bool uses_temporal_sampling;
    // The time variables (for timing the simulation in real time)
    time_t start, sim_start, sim_end, now, sim_finish, out_finish;
    time_t time_taken;
    // Active lineages stored as a row of datapoints
    vector<DataPoint> active;
    // Stores the point of the end of the active vector. 0 is reserved as null
    unsigned long endactive;
    // the maximum size of endactive
    unsigned long startendactive;
    // the maximum simulated number of individuals in the present day.
    unsigned long maxsimsize;
    // for create the link to the speciationcounter object which handles everything.
    Community community;
    // This might need to be updated for simulations that have large changes in maximum population size over time.
    // number of simulation num_steps
    long steps;
    // Maximum time to run for (in seconds)
    unsigned long maxtime;
    // number of generations passed,
    double generation;
    // The number of individuals per cell
    double deme;
    // The proportion of individuals to sample
    double deme_sample;
    // the speciation rate
    long double spec;
    // Path to output directory
    string out_directory;
    // sqlite3 object that stores all the data
    shared_ptr<SQLiteHandler> database;
    // only set to true if the simulation has finished, otherwise will be false.
    bool sim_complete;
    // set to true when variables are imported
    bool has_imported_vars;
// If sql database is written first to memory, then need another object to contain the in-memory database.
#ifdef sql_ram
    SQLiteHandler outdatabase;
#endif
    // Create the step object that will be retained for the whole simulation.
    // Does not need saving on simulation pause.
    Step this_step;
    string sql_output_database;
    // If true, means the command-line imports were under the (deprecated) fullmode.
    bool bFullMode;
    // If true, the simulation is to be resumed.
    bool bResume;
    // If true, a config file contains the simulation variables.
    bool bConfig;
    // If true, simulation can be resumed.
    bool has_paused, has_imported_pause;
    // Should always be false in the base class
    bool bIsProtracted;
    // variable for storing the paused sim location if files have been moved during paused/resumed simulations!
    string pause_sim_directory;
public:
    Tree() : data(make_shared<vector<TreeNode>>()), enddata(0), sim_parameters(make_shared<SimParameters>()),
             NR(make_shared<RNGController>()), speciation_rates(), seeded(false),
             seed(-1), job_type(-1), times_file("null"), reference_times(), uses_temporal_sampling(false),
             start(0), sim_start(0), sim_end(0), now(0), sim_finish(0), out_finish(0), time_taken(0), active(),
             endactive(0), startendactive(0), maxsimsize(0), community(data), steps(0), maxtime(0), generation(0.0),
             deme(0.0), deme_sample(0.0), spec(0.0), out_directory(""), database(make_shared<SQLiteHandler>()),
             sim_complete(false),
             has_imported_vars(false),
#ifdef sql_ram
             outdatabase(),
#endif //sql_ram
             this_step(), sql_output_database("null"), bFullMode(false), bResume(false), bConfig(true),
             has_paused(false), has_imported_pause(false), bIsProtracted(false), pause_sim_directory("null")
    {

    }

    virtual ~Tree()
    {
        database->close();
#ifdef sql_ram
        outdatabase.close();
#endif
    }

    /**
     * @brief Import the simulation variables from the command line structure.
     *
     * This function parses the simulation variables, imports them from the config file,
     * checks that the input files exist and checks for any paused simulations. The flags are then set correctly,
     * meaning that setup() and runSimulation() can be run immediately afterwards.

     * @param configfile the path to the config file containing parameters to parse.
     */
    void importSimulationVariables(const string &configfile);

    /**
     * @brief Import the simulation variables from a ConfigOption.
     *
     * This function parses the simulation variables, imports them (from either the command line or a config file),
     * checks that the input files exist and checks for any paused simulations. The flags are then set correctly,
     * meaning that setup() and runSim() can be run immediately afterwards.

     * @param config the set of config parameters to import
     */
    void importSimulationVariables(ConfigParser config);

    /**
     * @brief Runs the basic file existence checks.
     * Checks for paused simulations and file existence.
     */
    virtual void runFileChecks();

    /**
     * @brief Resets all the simulation variables.
     */
    void wipeSimulationVariables();

    /**
     * @brief Sets up the simulation parameters from the one provided.
     *
     * Intended for usage with metacommunity application. No output directory is expected.
     * @param sim_parameters_in the simulation parameters to set up the simulation with
     */
    void internalSetup(shared_ptr<SimParameters> sim_parameters_in);

    /**
     * @brief Asserts that the output directory is not null and exists. If it doesn't exist, it attempts to create it.
     * @throws Fatal_Exception if the output directory creation fails
     * @return true if output creates successfully
     */
    bool checkOutputDirectory();

    /**
     * @brief Checks for existing paused simulations to resume from
     * Sets bPaused if there are.
     *
     * This version uses the default values read from the config file.
     */
    void checkSims();

    /**
     * @brief Checks for existing paused simulations to resume from
     *
     * Sets bPaused if there are.
     *
     * This version uses the values supplied to the function directly
     *
     * @param output_dir the output directory to check for
     * @param seed the seed for paused sims
     * @param job_type the job_type for paused sims
     */
    void checkSims(string output_dir, long seed, long job_type);

    /**
     * @brief Move the parameters from the sim_parameters object to their relevant parameters.
     */
    virtual void setParameters();

    /**
     * @brief Sets the protracted variables
     * @param speciation_gen_min the minimum number of generations to have passed before speciation is allowed
     * @param speciation_gen_max the maximum number of generations a lineage can exist for before it is speciated.
     */
    virtual void setProtractedVariables(double speciation_gen_min, double speciation_gen_max);

    /**
     * @brief Gets the has_paused variable for resuming sims.
     * @return if the simulation has paused
     */
    bool hasPaused();

    /**
     * @brief Gets the map autocorrel times
     */
    vector<double> getTemporalSampling();

    /**
    * @brief Getter for the simulation seed.
    * @return Returns the seed
    */
    virtual long long getSeed();

    /**
    * @brief Gets the job type for the simulation.
    * This is a reference number for the jobs.
    * @return Returns the job type
    */
    virtual long long getJobType();

    /**
     * @brief Sets the simulation seed for the random number generator.
     *
     * This function should only be called once.
     *
     * The seed is set within the NR  object. This will be fixed for the simulation and is only performed once.
     *
     * @param seed_in the desired seed to set for the simulation
     */
    void setSeed(long long seed_in);

    /**
     * @brief Gets the initial number of individuals
     * @return the number of initial individuals to simulate
     */
    virtual unsigned long getInitialCount();

    /**
     * @brief Sets the sizes of grid, active and data, based on the number of individuals counted from the samplemask.
     * @return a count of the number of individuals that exist in the simulation
     */
    unsigned long setObjectSizes();

    /**
     * @brief The setup function for generating the simulation objects.
     *
     * The simulation parameters are set from comargs using setParameters().
     * Generates and fills the active and grid objects as well as importing all the maps from the supplied files.
     */
    virtual void setup();

    /**
     * @brief Sets the starting values for required parameters.
     */
    void setInitialValues();

    /**
     * @brief Sets the variables at the start of a simulation for temporary data.
     *
     * This is not the main set-up routine, which creates the permanent data structures including maps,
     * the coalescence tree and active lineage listings.
     */
    void setSimStartVariables();

    /**
     * @brief Prints the statement for the setup initiation.
     *
     * This is stored in a separate function so that it can be called in isolation by child classes.
     */
    void printSetup();

    /**
     * @brief Sets the temporal sampling points from the time config file.
     */
    void setTimes();

    /**
     * @brief Determines the speciation rates to apply and then applies them to the coalescence tree post-simulation.
     *
     * Detects speciation rates from the config files supplied.
     */
    void determineSpeciationRates();

    /**
     * @brief Adds the speciation rates to those to be applied.
     */
    void addSpeciationRates(vector<long double> spec_rates_in);

    /**
     * @brief Assigns the objects sizes in memory and fills with the starting lineages.
     */
    void generateObjects();

    /**
     * @brief Fills the active and data objects with the starting lineages.
     * @param initial_count the number of individuals expected to exist
     */
    virtual unsigned long fillObjects(const unsigned long &initial_count);

    /**
     * @brief Run the entire simulation given the start conditions already defined by setup()
     *
     * Setup is assumed to have been run already. This function is the main function containing the main loop of the
     * simulation.
     * At the end of the simulation, returns true if the simulation is complete, false otherwise.
     */
    virtual bool runSimulation();

    /**
     * @brief Writes to the console that the simulation is beginning
     */
    void writeSimStartToConsole();

    /**
     * @brief Write the step counter to console. This function should only be called in debugging mode.
     */
    void writeStepToConsole();

    /**
     * @brief Increments the generation counter and step references
     */
    virtual void incrementGeneration();

    /**
     * @brief Chooses a random lineage from active.
     *
     * The index of the random lineage is stored in this_step, as chosen.
     * Also records the required variables for the step process, like x, y position.
     */
    void chooseRandomLineage();

    /**
     * @brief Updates the coalescence variables in the step object.
     */
    virtual void updateStepCoalescenceVariables();

    /**
     * @brief Speciation to supplied lineage.
     *
     * Also calls the removeOldPos() and switchPositions() functions for removing the lineage out of active reference.
     * @param chosen
     */
    void speciation(const unsigned long &chosen);

    /**
     * @brief Performs the actual speciation.
     * @param data_position the position in the array of TreeNodes for this lineage
     */
    virtual void speciateLineage(const unsigned long &data_position);

    /**
     * @brief Removes the old position within active.
     * @param chosen the desired active reference to remove from the grid.
     */
    virtual void removeOldPosition(const unsigned long &chosen);

    /**
     * @brief Switches the chosen position with the endactive position.
     * @param chosen the chosen lineage to switch with endactive.
     */
    virtual void switchPositions(const unsigned long &chosen);

    /**
     * @brief Calculates the next step for the simulation.
     */
    virtual void calcNextStep();

    /**
     * @brief Calculates the speciation probability from the random number, speciation rate and number of generations a
     * lineage has existed for.
     * @param random_number the generated random number from 0-1
     * @param speciation_rate the speciation rate to be applied
     * @param no_generations the number of generations a lineage has existed for
     * @return if true, speciation has occured
     */
    virtual bool calcSpeciation(const long double &random_number,
                                const long double &speciation_rate,
                                const unsigned long &no_generations);

    /**
    * @brief  Perform the coalescence between lineages. Once coalesced, lineages are removed from the active scope.
    * @param chosen the chosen lineage for coalescence
    * @param coalchosen the target lineage for coalscence
    */
    void coalescenceEvent(const unsigned long &chosen, unsigned long &coalchosen);

    /**
     * @brief Checks if the number of lineages should be expanded at another sample point
     */
    void checkTimeUpdate();

    /**
     * @brief Adds the required lineages at the generation time.
     * @param generation_in the generation to add lineages at
     */
    virtual void addLineages(double generation_in);

    /**
     * @brief Randomly checks if a lineage should be added, based on the proportion added.
     * @param proportion_added the proportion of lineages that should be added
     * @return true if the lineage should be added, false otherwise
     */
    bool checkProportionAdded(const double &proportion_added);

    /**
     * @brief Checks the size of the main active and data objects is large enough
     * @param req_data the required data object size
     * @param req_active the required active object size
     */
    void checkSimSize(unsigned long req_data, unsigned long req_active);

    /**
     * @brief Sets the active reference to a tip, if it isn't one already. Otherwise, creates a new tip for the new
     * generation time.
     * @param tmp_active the reference in active
     * @param generation_in the generation to set for the new lineage
     * @param data_added vector of lineages to add to data
     */
    void makeTip(const unsigned long &tmp_active, const double &generation_in, vector<TreeNode> &data_added);

    /**
     * @brief Creates a new reference in data containing the tip with a new generation counter.
     * @param i the reference in active of the lineage to make a tip
     * @param generationin the generation to make the lineage a tip at
     * @param data_added vector of lineages to add to the data object
     */
    void convertTip(unsigned long i, double generationin, vector<TreeNode> &data_added);

    /**
     * @brief Finalises the simulation, and performs the correct tasks depending if the sim has been paused or finished.
     * @return
     */
    bool stopSimulation();

    /**
     * @brief Applies the given speciation rate to the tree.
     * @param sr the required speciation rate
     * @param t the required time of speciation
     */
    void applySpecRate(long double sr, double t);

    /**
     * @brief Applies the given speciation rate to the tree, but does not output to a file. Instead returns a pointer
     * to the nodes object.
     * @param sr the required speciation rate
     * @param t the required time of speciation
     */
    void applySpecRateInternal(long double sr, double t);

    /**
     * @brief Gets the sorted cumulative species abundances from the contained TreeList.
     *
     * For use with metacommunity applications
     * @return row of cumulative species abundances
     */
    shared_ptr<vector<unsigned long>> getCumulativeAbundances();

    /**
     * @brief Gets the species abundances from the internal tree.
     * @param community_reference the community reference
     * @return the species abundances
     */
    shared_ptr<map<unsigned long, unsigned long>> getSpeciesAbundances(const unsigned long &community_reference);

    /**
     * @brief Gets the species abundances from the internal tree.
     * @return the species abundances
     */
    shared_ptr<vector<unsigned long>> getSpeciesAbundances();

    /**
     * @brief Sets up Community member to point to the same output database as the simulation.
     * @return the protracted speciation parameters to set up
     */
    ProtractedSpeciationParameters setupCommunity();

    /**
     * @brief Sets up the generation of the tree object.
     * @param sr the required speciation rate
     * @param t the required time of speciation
     */
    void setupCommunityCalculation(long double sr, double t);

    /**
     * @brief Overloaded version of applySpecRates for the default generation (0.0).
     * @param sr the speciation rate to apply to the tree
     */
    void applySpecRate(long double sr);

    /**
     * @brief Applies multiple speciation rates to the coalescence tree, ignoring repeated speciation rates.
     *
     * Speciation rates are read from the speciation_rates object, which should have already been calculated.
     */
    void applyMultipleRates();

    /**
     * @brief Gets the protractedness of the simulation.
     * Overridden by protracted child classes.
     * @return
     */
    virtual bool getProtracted();

    /**
     * @brief Gets the protracted variables and returns them as a single, newline separated string.
     * This method is intended to be overridden in derived classes.
     * @return string containing the protracted variables, separated by newlines.
     */
    virtual string getProtractedVariables();

    /**
     * @brief Gets the minimum number of generations a lineage must exist.
     *
     * Without overriding this function, should always return 0.0.
     * @return double the number of generations a lineage must exist
     */
    virtual double getProtractedGenerationMin();

    /**
     * @brief Gets the maximum number of generations a lineage can exist.
     *
     * Without overriding this function, should always return 0.0 (no maximum).
     *
     * @return double the number of generations a lineage must exist
     */
    virtual double getProtractedGenerationMax();

    /**
     * @brief Copy the in-memory database to file.
     *
     * This function should not be called if the database is already opened on disc, and won't do anything if it is.
     */
    void sqlOutput();

    /**
     * @brief Creates the output data objects and outputs important simulation data to a database file.
     */
    void createAndOutputData();

    /**
     * @brief Outputs the simulation data (including the coalescence tree) to an sql database.
     */
    void outputData();

    /**
     * @brief Sort and process the species list so that the useful information can be extracted from it.
     */
    void sortData();

    /**
     * @brief Writes the times to the terminal for simulation information.
     */
    void writeTimes();

    /**
     * @brief Opens a connection to the in-memory database, or the on-disk database, depending on the compilation
     * options.
     */
    void openSQLDatabase();

    /**
     * @brief Generates the SQL database file from the full simulation data. This allows for greater analysis of the
     * data after completion of the simulation
     */
    void sqlCreate();

    /**
     * @brief Creates the output directory, if it doesn't already exist, and deletes any existing database with the
     * output name.
     */
    void setupOutputDirectory();

    /**
     * @brief Creates the SIMULATION_PARAMETERS table in the SQL database.
     */
    void sqlCreateSimulationParameters();

    /**
     * @brief Creates a string containing the SQL insertion statement for the simulation parameters.
     * @return string containing the SQL insertion statement
     */
    virtual string simulationParametersSqlInsertion();

    /**
     * @brief Outputs the protracted variables to a string.
     *
     * This function is intended to be overridden by derived classes.
     * It is intended the output is used for writing to SQL databases.
     *
     * @return string containing a list of the protracted speciation variables.
     */
    virtual string protractedVarsToString();

    /**
     * @brief Pause the simulation and dump data from memory.
     */
    virtual void simPause();

    /**
     * @brief Checks the output folder exists and initiates the pause.
     * @return the output file stream to save objects to
     */
    shared_ptr<ofstream> initiatePause();

    /**
     * @brief Saves the main simulation variables to file.
     * @param out the output file stream to save the object to
     */
    void dumpMain(shared_ptr<ofstream> out);

    /**
     * @brief Saves the active object to file.
     * @param out the output file stream to save the object to
     */
    void dumpActive(shared_ptr<ofstream> out);

    /**
     * @brief Saves the data object to file.
     * @param out the output file stream to save the object to
     */
    void dumpData(shared_ptr<ofstream> out);

    /**
     * @brief Completes the pause routine and outputs the sql dump.
     * @param out the output stream to close up
     */
    void completePause(shared_ptr<ofstream> out);

    /**
     * @brief Sets the resume variables so that the simulation can be resumed.
     *
     * The pause directory can be the same as the output directory if it is desirable to save to the same location.

     * @param pausedir the directory containing the pause folder for resuming the simulation
     * @param outdir the directory to write simulation output to
     * @param seed the simulation seed
     * @param job_type the simulation job reference number
     * @param new_max_time the maximum simulation time to run for in seconds (0 keeps old simulation max time)
     */
    void setResumeParameters(string pausedir, string outdir, unsigned long seed, unsigned long job_type,
                             unsigned long new_max_time);

    /**
     * @brief Sets the resume variables to the defaults
     */
    void setResumeParameters();

    shared_ptr<ifstream> openSaveFile();

    /**
     * @brief Loads the main simulation parameters from the save file into memory.
     */
    virtual void loadMainSave(shared_ptr<ifstream> in1);

    /**
     * @brief Loads the data object from the save file into memory.
     */
    void loadDataSave(shared_ptr<ifstream> in1);

    /**
     * @brief Loads the active object from the save file into memory.
     */
    void loadActiveSave(shared_ptr<ifstream> in1);

    /**
     * @brief Checks for resuming and prints to the terminal
     */
    void initiateResume();

    /**
     * @brief Resumes the simulation from a previous state.
     *
     * Reads in the parameters and objects from file and re-starts the simulation.
     */
    virtual void simResume();

#ifdef DEBUG

    /**
     * @brief Validates all lineages have been set up correctly. This may take considerable time for larger simulations.
     *
     * This function is only relevant in debug mode.
     */
    virtual void validateLineages();

    /**
     * @brief Runs the debug checks at the end of each step.
     *
     * Should only be run if compiled with the DEBUG flag.
     */
    virtual void debugEndStep();

    /**
     * @brief Runs the debug checks upon a coalescence event.
     *
     * Should only be run if compiled with the DEBUG flag.
     */
    void debugCoalescence();

    /**
     * @brief Run checks at the end of each cycle which make certain the step has been successful.
     * @param chosen the chosen lineage to check
     * @param coalchosen the lineage which is coalescing with the chosen lineage which we are also required to check
     */
    virtual void runChecks(const unsigned long &chosen, const unsigned long &coalchosen);

    /**
     * @brief Runs a mini check on the chosen lineage.
     * @param chosen the lineage to check for
     */
    void miniCheck(const unsigned long &chosen);
#endif // DEBUG
};

#endif //TREE_H
