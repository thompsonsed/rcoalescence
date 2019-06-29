//This file is part of necsim project which is released under MIT license.
//See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.
/**
 * @file Community.h
 *
 * @brief Contains various objects used for reconstructing the coalescence tree after simulations are complete.
 *
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */

#ifndef TREELIST
#define TREELIST

#include <cmath>
#include <sqlite3.h>
#include <cstring>
#include <stdexcept>
#include <string>
# include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <set>
#include <utility>
#include <memory>

#ifdef WIN_INSTALL
#include <windows.h>
#define sleep Sleep
#endif

#include "TreeNode.h"
#include "Matrix.h"
#include "DataMask.h"
#include "parameters.h"
#include "SpecSimParameters.h"
#include "SQLiteHandler.h"

using namespace std;
using std::string;

/**
 * @brief Checks whether speciation has occured for the provided parameters. 
 * Provided here for ease of use when bug-fixing.
 * @param random_number the random number associated with a lineage
 * @param speciation_rate the global speciation rate
 * @param number_of_generations the number of generations the lineage has existed
 * @return bool the speciation state of the lineage
 */
bool checkSpeciation(const long double &random_number, const long double &speciation_rate,
                     const unsigned long &no_generations);

/**
 * @brief Contains the information needed for defining a fragment.
 */
struct Fragment
{
    // the name for the fragment (for reference purposes)
    string name;
    // coordinates for the extremes of the site
    long x_east, x_west, y_north, y_south;
    // the number of lineages in the fragment.
    unsigned long num;
    double area;
};

/**
 * @brief A child of the Matrix class as booleans.
 * Used for determining where to sample species from.
 */
class Samplematrix : public DataMask
{
private:
    bool bIsNull;
    bool bIsFragment;
    Fragment fragment;
public:
    /**
     * @brief Inherit construction from the Matrix class, but also set the booleans.
     */
    Samplematrix();

    //	/**
    //	 * @brief Returns the value at the x,y position.
    //	 * This is used for testing purposes only.
    //	 * @param xval the x coordinate.
    //	 * @param yval the y coordinate
    //	 * @param xwrap the x wrapping
    //	 * @param ywrap the y wrapping
    //	 * @return the value at x,y.
    //	 */
    bool getTestVal(unsigned long xval, unsigned long yval, long xwrap, long ywrap);

    /**
     * @brief Returns the value at the x,y position, with the given x and y wrap.
     * Also checks whether or not the map is set to null, or whether the value comes from within a fragment.
     * @param x1 the x coordinate.
     * @param y1 the y coordinate
     * @param xwrap the x wrapping
//	 * @param ywrap the y wrapping
     * @return the value at x,y.
     */
    bool getMaskVal(unsigned long x1, unsigned long y1, long x_wrap, long y_wrap);

    /**
     * @brief Set the fragment for the samplemask to some calculated fragment.
     * This can be set multiple times.
     * @param fragment_in the Fragment to set the samplemask to.
     */
    void setFragment(Fragment &fragment_in);

    /**
     * @brief Removes the fragment.
     */
    void removeFragment();
};

/**
 * @brief A class to contain the tree object lineages and reconstructing the coalescence tree.
 * Contains functions for calculating the number of species for a given speciation rate,
 * outputting spatial data and generating species abundance distributions.
 * Requires a link to the SQLite database from simulation output,
 * and produces results within the same database file.
 */
class Community
{
protected:
    bool in_mem; // boolean for whether the database is in memory or not.
    bool database_set; // boolean for whether the database has been set already.
    shared_ptr<SQLiteHandler> database; // stores the in-memory database connection.
    bool bSqlConnection; // true if the data connection has been established.
    shared_ptr<vector<TreeNode>> nodes; // in older versions this was called species_id_list.
    shared_ptr<vector<unsigned long>> species_abundances;
    unsigned long iSpecies;
    bool has_imported_samplemask; // checks whether the samplemask has already been imported.
    bool has_imported_data; // checks whether the main sim data has been imported.
    Samplematrix samplemask; // the samplemask object for defining the areas we want to sample from.
    vector<Fragment> fragments; // a vector of fragments for storing each fragment's coordinates.
    shared_ptr<CommunityParameters> current_community_parameters;
    shared_ptr<MetacommunityParameters> current_metacommunity_parameters;
    // the minimum speciation rate the original simulation was run with
    // this is read from the database SIMULATION_PARAMETERS table
    long double min_spec_rate;
    // The dimensions of the sample grid size.
    unsigned long grid_x_size, grid_y_size;
    // The dimensions of the original sample map file
    unsigned long samplemask_x_size, samplemask_y_size, samplemask_x_offset, samplemask_y_offset;
    // Vector containing past speciation rates
    CommunitiesArray past_communities;
    MetacommunitiesArray past_metacommunities;
    // Protracted speciation current_metacommunity_parameters
    bool protracted;
    ProtractedSpeciationParameters minimum_protracted_parameters;
    ProtractedSpeciationParameters applied_protracted_parameters;
    unsigned long max_species_id, max_fragment_id, max_locations_id;
    // Does not need to be stored during simulation pause
    shared_ptr<SpecSimParameters> spec_sim_parameters;
public:

    /**
     * @brief Contructor for the community linking to Treenode list.
     * @param r Row of TreeNode objects to link to.
     */
    explicit Community(shared_ptr<vector<TreeNode>> r) : in_mem(false), database_set(false),
                                                         database(make_shared<SQLiteHandler>()), bSqlConnection(false),
                                                         nodes(std::move(r)),
                                                         species_abundances(make_shared<vector<unsigned long>>()),
                                                         iSpecies(0), has_imported_samplemask(false),
                                                         has_imported_data(false), samplemask(), fragments(),
                                                         current_community_parameters(
                                                                 make_shared<CommunityParameters>()),
                                                         current_metacommunity_parameters(
                                                                 make_shared<MetacommunityParameters>()),
                                                         min_spec_rate(0.0),
                                                         grid_x_size(0), grid_y_size(0), samplemask_x_size(0),
                                                         samplemask_y_size(0), samplemask_x_offset(0),
                                                         samplemask_y_offset(0), past_communities(),
                                                         past_metacommunities(), protracted(false),
                                                         minimum_protracted_parameters(),
                                                         applied_protracted_parameters(), max_species_id(0),
                                                         max_fragment_id(0), max_locations_id(0),
                                                         spec_sim_parameters(make_shared<SpecSimParameters>())
    {

    }

    /**
     * @brief Default constructor
     */
    Community() : Community(make_shared<vector<TreeNode>>())
    {
    }

    /**
    * @brief Default destructor
    */
    virtual ~Community()
    {
        if(nodes)
        {
            nodes.reset();
        }
        closeSqlConnection();
    }

    /**
     * @brief Set the nodes object to the input Row of Treenode objects.
     * @param l the Row of Treenode objects to link to.
     */
    void setList(shared_ptr<vector<TreeNode>> l);

    /**
     * @brief Sets up the community from a set of simulation parameters and the sqlite3 database connection.
     * @param sim_parameters parameters to use for setting up the community
     * @param database points to the database object to open
     * @return
     */
    ProtractedSpeciationParameters setupInternal(shared_ptr<SimParameters> sim_parameters,
                                                 shared_ptr<SQLiteHandler> database);

    /**
     * @brief Sets the database object for the sqlite functions.
     * @param dbin the sqlite3 input database.
     */
    void setDatabase(shared_ptr<SQLiteHandler> dbin);

    /**
     * @brief Get the boolean of whether the data has been imported yet.
     * @return true if database has been imported.
     */
    bool hasImportedData();

    /**
     * @brief Get the minimum speciation rate the simulation was originally run with.
     * This value is read in from the SIMULATION_PARAMETERS table in the database file.
     * @return the minimum speciation rate.
     */
    long double getMinimumSpeciation();

    /**
     * @brief Imports the samplemask if it hasn't already been imported.
     * @param sSamplemask the path to the samplemask file.
     */
    void importSamplemask(string sSamplemask);

    /**
     * @brief Counts the number of species that have speciated currently on the tree.
     *
     * @return the number of species
     */
    unsigned long countSpecies();

    /**
     * @brief Calculate the number of species in the list for the parameters in the current_community_parameters object.
     * This is the main function which reconstructs the coalescence tree.
     * Each Treenode object will end having its existence value set correctly after a call to this function.

     * @return the number of species present.
     */
    unsigned long calculateCoalescencetree();

    /**
     * @brief Speciates TreeNode and updates the species count.
     *
     * For systems which are not using a metacommunity, this function will just perform basic speciation.
     *
     * @note species_list is not updated in unless the function is overridden for metacommunity application.
     * @param species_count the total number of species currently in the community
     * @param treenode pointer to the TreeNode object for this lineage
     * @param species_list the set of all species ids.
     */
    virtual void addSpecies(unsigned long &species_count, TreeNode* treenode, set<unsigned long> &species_list);

    /**
     * @brief Calculates the species abundance of the dataset.
     * The species abundances will be with rOut after a call do this function.
     * If a samplemask has been applied, only lineages which originally existed in the samplemask will be counted.
     */
    void calcSpeciesAbundance();

    /**
     * @brief Resets the entire tree.
     * Sets existance to false, speciation to false and removes any species ID.
      */
    virtual void resetTree();

    /**
     * @brief Opens the connection to the sql database file
     * Note that this imports the database to memory, so functionality should be changed for extremely large database
     * files.
     * @param input_file the sql database output from a necsim simulation.
     */
    void openSqlConnection(string input_file);

    /**
     * @brief Safely destroys the SQL connection.
     */
    void closeSqlConnection();

    /**
     * @brief Opens a connection to an in-memory database. This will eventually be written to the output file.
     */
    void setInternalDatabase();

    /**
     * @brief Internally sets the file referencing, data import and sql connection flags to true, for allowing checks
     * to pass from internal object creation (so no external files are needed)
     */
    void internalOption();

    /**
     * @brief Imports the data from the desired SQL database object into the array.
     * @note Opens the sql connection if it has not already been opened.
     * @note If nodes is not of length 0, this function does nothing. This is so that any in-memory data is not
     * overwritten.
     * @param inputfile the path to the input SQLite database.
     */
    void importData(string inputfile);

    /**
     * @brief Sets the simulation parameters from a SimParameters object.
     * @param sim_parameters pointer to the SimParameters object to set from
     */
    void setSimParameters(shared_ptr<SimParameters> sim_parameters);

    /**
     * @brief Sets the speciation parameters from a SpecSimParameters object.
     * @note Be careful of the difference between this function and setSimParameters (which sets the main simulation
     * parameters).
     * @param spec_sim_parameters pointer to the SpecSimParameters object containing community speciation parameters
     */
    void setSpecSimParameters(shared_ptr<SpecSimParameters> spec_sim_parameters);

    /**
     * @brief Imports the simulation parameters by reading the SIMULATION_PARAMETERS table in the provided file.
     * This imports the grid_x_size, grid_y_size (which should also be the sample map dimensions) and the minimum
     * speciation rate.
     *
     * @note Opens the sql connection if it has not already been opened.
     *
     * @note If bDataImport has already been set, no operation is performed.
     *
     * @param file the sqlite database simulation output which will be used for coalescence tree generation.
     */
    void importSimParameters(string file);

    /**
     * @brief Forces the sim_complete parameter to be true in SIMULATION_PARAMETERS.
     * Used when speciating all remaining lineages to force a simulation to completion.
     */
    void forceSimCompleteParameter();

    /**
     * @brief Gets if the database has been set yet.
     * @return true if the database is already set
     */
    bool isSetDatabase();

    /**
     * @brief Gets the maximum species abundance ID from the database and stores it in the max_species_id variable.
     * @note Does not check for SPECIES_ABUNDANCES existence and will throw an error if it cannot access it
     */
    void getMaxSpeciesAbundancesID();

    /**
     * @brief Changes the rOut object so that values represent cummulative species abundances.
     *
     * Allows binary sort on rOut (much faster) and the previous rOut value can be obtained by
     * value = rOut[i] - rOut[i-1]
     * @return pointer to sorted Row of species abundances
     */
    shared_ptr<vector<unsigned long>> getCumulativeAbundances();

    /**
     * @brief Returns the row_out object, which should contain species abundances or cumulative abundances
     * @note Does not recalculate species abundances, so if getCumulativeAbundances has been called, will return the
     * cumulative species abundances instead.
     * @note Returns a copy, so could cause problems for extremely large simulations with immense numbers of species.
     * @return row_out, the species abundances, or the cumulative abundances if getCumulativeAbundances has been called
     */
    shared_ptr<vector<unsigned long>> getRowOut();

    /**
     * @brief Gets the number of species in the most recent calculation.
     * @return the number of species in the most recent calculation
     */
    unsigned long getSpeciesNumber();

    /**
     * @brief Gets the maximum fragment abundance ID from the database and stores it in the max_fragment_id variable.
     * @note Does not check for FRAGMENT_ABUNDANCES existence and will throw an error if it cannot access it
     */
    void getMaxFragmentAbundancesID();

    /**
     * @brief Gets the maximum species locations ID from the database and stores it in the max_locations_id variable.
     * @note Does not check for SPECIES_LOCATIONS existence and will throw an error if it cannot access it

     */
    void getMaxSpeciesLocationsID();

    /**
     * @brief Sets the protracted parameters for application of protracted speciation
     *
     * This overloaded version is for setting protracted parameters before a full simulation has been outputted (i.e.
     * immediately after completion of the simulation).
     *
     * @param protracted_params protracted speciation parameters to add
     */
    void setProtractedParameters(const ProtractedSpeciationParameters &protracted_params);

    /**
     * @brief Overrides the protracted parameters for the Community object
     * @param protracted_params the protracted parameters to override with
     */
    void overrideProtractedParameters(const ProtractedSpeciationParameters &protracted_params);

    /**
     * @brief Sets the protracted boolean to the input.
     *
     * @param protracted_in the protracted boolean to set
     */
    void setProtracted(bool protracted_in);

    /**
     * @brief Creates a new table in the database file and outputs the database object to the same file as the input
     * file. Calculates the community structure for the set of community parameters in current_community_parameters.
     *
     * The new SPECIES_ABUNDANCES table contains the species abundance distribution for the whole samplemask.
     * A similar tabe FRAGMENT_ABUNDANCES is generated by createFragmentDatabase() if specified via the command line
     * parameters.
     */
    void createDatabase();

    /**
     * @brief Calculates the coalescence tree and generates species abundances.
     */
    void generateBiodiversity();

    /**
     * @brief Outputs the species abundances into the database
     */
    void outputSpeciesAbundances();

    /**
     * @brief Checks if calculations with the given set of parameters has already been performed.
     * @param speciation_rate the speciation rate to check for
     * @param time the time to check for
     * @param fragments if true, checks fragments have been used
     * @param metacommunity_size the metacommunity size to check for
     * @param metacommunity_speciation_rate the metacommunity speciation rate to check for
     * @param proc_parameters protracted speciation parameters to add
     * @return
     */
    bool checkCalculationsPerformed(const long double &speciation_rate, const double &time, const bool &fragments,
                                    const MetacommunityParameters &metacomm_parameters,
                                    const ProtractedSpeciationParameters &proc_parameters);

    /**
     * @brief Adds a performed calculation to the lists of calculations. Also sets the current_community_parameters
     * pointer to the set of parameters to be applied.
     *
     * @param speciation_rate the speciation rate of the performed calculation
     * @param time the time of the performed calculation
     * @param fragments if true, fragments were used
     * @param metacommunity_size the metacommunity size of the performed calculation
     * @param metacommunity_speciation_rate the metacommunity speciation rate of the performed calculation
     * @param protracted_parameters protracted speciation parameters to add
     */
    void addCalculationPerformed(const long double &speciation_rate, const double &time, const bool &fragments,
                                 const MetacommunityParameters &metacomm_parameters,
                                 const ProtractedSpeciationParameters &protracted_parameters);

    /**
     * @brief Creates a new table in the database file and outputs the database object to the same file as the input file.
     * Essentially creates a species abundance distribution (as in createDatabase()), but for the specified fragment
     * within the samplemask.
     * @param f the Fragment to sample from.
     */
    void createFragmentDatabase(const Fragment &f);

    /**
     * @brief Output the database from memory to the database file.
     * Most of the time, it is desirable for the outputfile to be the same path
     * as the input file and will write to the same object.
     */
    void exportDatabase();

    /**
     * @brief Checks for the current CommunityParameters reference in the SPECIES_LOCATIONS table.
     * @return true if the reference exists in the SPECIES_LOCATIONS table
     */
    bool checkSpeciesLocationsReference();

    /**
     * @brief Checks for the current CommunityParameters reference in the SPECIES_ABUNDANCES table.
     * @return return true if the reference exists in the SPECIES_LOCATIONS table
     */
    bool checkSpeciesAbundancesReference();

    /**
     * @brief Record the full spatial data.
     * Creates a new table, SPECIES_LOCATIONS containing every species and their parameters.
     * This allows for more in-depth analysis to be performed if necessary.
     */
    void recordSpatial();

    /**
     * @brief Calculates the limits of each fragment in the sample map and adds it to the vector of fragments.
     * If the fragment_file is null, then the program will attempt to calculate fragments from the map.
     * @bug Only rectangular fragments will be detected. Problems will also be encountered for adjacent fragments.
     * @param fragment_file the fragment file to read from.
     */
    void calcFragments(string fragment_file);

    /**
     * @brief Calculate species abundances for each fragment, and call createFragmentDatabase() for each Fragment.
     */
    void applyFragments();

    /**
     * @brief Gets the previous calculations that have already been performed.
     */
    void getPreviousCalcs();

    /**
     * @brief Gets the unique community references from the SQL database
     * @return a vector containing the unique references
     */
    vector<unsigned long> getUniqueCommunityRefs();

    /**
     * @brief Gets the unique metacommunity reference from the SQL database
     * @return a vector containing the unique references
     */
    vector<unsigned long> getUniqueMetacommunityRefs();

    /**
     * @brief Write all performed calculations to the output database.
     */
    void writeNewCommunityParameters();

    /**
     * @brief Write all performed calculations to the output database
     */
    void writeNewMetacommunityParameters();

    /**
     * @brief Creates a new table, SPECIES_LIST in the output database.
     */
    void createSpeciesList();

    /**
     * @brief Drops the SPECIES_LIST table from the database.
     */
    void deleteSpeciesList();

    /**
     * @brief Writes the coalescence tree to a table called SPECIES_LIST.
     */
    void writeSpeciesList(const unsigned long &enddata);

    /**
     * @brief Updates the fragments tag on those simulations which now have had fragments added.
     */
    void updateCommunityParameters();

    /**
     * @brief Prints speciation rates to terminal
     */
    void writeSpeciationRates();

    /**
     * @brief Calculates the coalescence tree for each set of parameters in speciation_parameters;
     */
    void calculateTree();

    /**
     * @brief Outputs the data to the SQL database.
     */
    void output();

    /**
     * @brief Prints the application times.
     * @param tStart the start time
     * @param tEnd the end time
     */
    void printEndTimes(time_t tStart, time_t tEnd);

    /**
     * @brief Apply the given speciation parameters to the coalescence tree.
     * Overridden for metacommunity application.
     * @param sp speciation parameters to apply, including speciation rate, times and spatial sampling procedure
     */
    void apply(shared_ptr<SpecSimParameters> sp);

    /**
     * @brief Applies the given speciation parameters to the coalescence tree, but does not write the output.
     * @param sp speciation parameters to apply, including speciation rate, times and spatial sampling procedure
     */
    void applyNoOutput(shared_ptr<SpecSimParameters> sp);

    /**
     * @brief Applies the given speciation parameters to the coalescence tree, but does not write the output.
     * @param sp speciation parameters to apply, including speciation rate, times and spatial sampling procedure
     * @param tree_data the coalescence tree containing simulation data
     */
    virtual void applyNoOutput(shared_ptr<SpecSimParameters> sp, shared_ptr<vector<TreeNode>> tree_data);

    /**
     * @brief Sets up the community application by reading parameters and data.
     * @param sp the speciation parameters to use for generating the community
     * @param data the list of all nodes on the coalescence tree
     */
    void setupApplication(shared_ptr<SpecSimParameters> sp, shared_ptr<vector<TreeNode>> data)
    {
        spec_sim_parameters = sp;
        writeSpeciationRates();
        // Set up the objects
        setList(std::move(data));
        importSimParameters(sp->filename);
        importSamplemask(sp->samplemask);
        importData(sp->filename);
        getPreviousCalcs();
        if(sp->use_fragments)
        {
            calcFragments(sp->fragment_config_file);
            stringstream os;
            os << "Total fragments: " << fragments.size() << endl;
            writeInfo(os.str());
        }
        if(spec_sim_parameters->metacommunity_parameters.empty())
        {
            spec_sim_parameters->metacommunity_parameters.addNull();
            current_metacommunity_parameters = spec_sim_parameters->metacommunity_parameters.metacomm_parameters[0];
        }
        if(spec_sim_parameters->protracted_parameters.empty())
        {
            ProtractedSpeciationParameters tmp;
            spec_sim_parameters->protracted_parameters.emplace_back(tmp);
        }
    }

    /**
     * @brief Creates the coalescence tree for the given speciation parameters.
     * @param sp speciation parameters to apply, including speciation rate, times and spatial sampling procedure
     */
    void doApplication(shared_ptr<SpecSimParameters> sp);

    /**
     * @brief Creates the coalescence tree for the given speciation parameters.
     * @param sp speciation parameters to apply, including speciation rate, times and spatial sampling procedure
     * @param data the Row of TreeNodes that contains the coalescence tree.
     */
    void doApplication(shared_ptr<SpecSimParameters> sp, shared_ptr<vector<TreeNode>> data);

    /**
     * @brief Creates the coalescence tree for the given speciation parameters, using internal file referencing
     * to avoid any actual file creation.
     * @param sp speciation parameters to apply, including speciation rate, times and spatial sampling procedure
     * @param data the Row of TreeNodes that contains the coalescence tree.
     */
    void doApplicationInternal(shared_ptr<SpecSimParameters> sp, shared_ptr<vector<TreeNode>> data);

    /**
     * @brief Speciates the remaining lineages in an incomplete simulation to force it to appear complete.
     */
    void speciateRemainingLineages(const string &filename);

    /**
     * @brief Gets the species richness for the community reference from the database.
     * @param community_reference the community reference to obtain the species richness for
     * @return the number of species
     */
    unsigned long getSpeciesRichness(const unsigned long &community_reference);

    /**
     * @brief Gets the species abundances of the community
     * @param community_reference the reference of the desired community
     * @return a map of species ids to species abundances
     */
    shared_ptr<map<unsigned long, unsigned long>> getSpeciesAbundances(const unsigned long &community_reference);

    /**
     * @brief Gets the species abundance from the species_abundances internal object
     * @return the species abundances
     */
    shared_ptr<vector<unsigned long>> getSpeciesAbundances();

    /**
     * @brief Check if the database is a null pointer.
     * @return true if the database is a null pointer.
     */
    bool isDatabaseNullPtr();
};

#endif