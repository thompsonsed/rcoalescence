// This file is part of rcoalescence project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.
/**
 * @author Sam Thompson
 * @file RTree.h
 * @brief R wrapper for running non-spatial neutral models.
 *
 * Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */

#ifndef RCOALESCENCE_RTREE_H

#define RCOALESCENCE_RTREE_H

#include <string>

#ifndef CXX14_SUPPORT

#include "memory.h"

#else
#include <memory>
#endif

#ifndef WIN_INSTALL
#include <unistd.h>
#endif // WIN_INSTALL

#include <Rcpp.h>
//#include <Rcpp/r/headers.h>
//#include <Rcpp/DataFrame.h>
#include "necsim/SpecSimParameters.h"
#include "necsim/Tree.h"
#include "necsim/Metacommunity.h"
#include "RLogging.h"

using namespace std;
using namespace necsim;

namespace rcoalescence
{

    /**
     * @brief Class containing non-spatial non-protracted Tree simulations.
     */
    class RTreeSimulation : public virtual Tree
    {
    protected:
        shared_ptr<SpecSimParameters> spec_sim_parameters;
        Metacommunity metacommunity; // Provide both a Community and Metacommunity object for both use case
        bool multiple_output;
        bool has_outputted;
        bool has_run_setup;
        bool has_written_main_sim;
        bool uses_metacommunity;

    public:
        string output_database;

        RTreeSimulation();

        ~RTreeSimulation() override;
        
        /**
         * @brief Checks if the simulation parameters have already been set (i.e. the simulation has already been setup).
         * @return true if the simulation has been setup
         */
        bool checkHasSetup();
        
        /**
         * @brief Sets the main simulation parameters in the sim_parameters object.
         * @param task the job reference number, used for file referencing
         * @param seed_in the seed to set random number generation
         * @param output_directory_in the output directory
         * @param max_time_in the maximum time to simulate for
         * @param desired_specnum_in the desired number of species to aim towards (currently not functional)
         * @param times_list the file containing a list of temporal sampling points
         */
        void setKeyParameters(const long long &task,
                              const long long &seed_in,
                              const string &output_directory_in,
                              const unsigned long &max_time_in,
                              const unsigned long &desired_specnum_in,
                              vector<double> times_list);

        /**
         * @brief Sets the speciation parameters for the simulation in the sim_parameters object.
         * @param spec_in the minimum speciation rate to use
         */
        void setMinSpeciationRate(const long double &spec_in);

        /**
         * @brief Sets the protracted variables
         * @param protracted_in if true, simulates as a protracted simulation
         * @param min_speciation_gen_in the minimum speciation generation for protracted simulations
         * @param max_speciation_gen_in the maximum speciation generation for protracted simulations
         */
        void setSimulationProtractedParameters(const bool &protracted_in,
                                               const double &min_speciation_gen_in,
                                               const double &max_speciation_gen_in);

        /**
         * @brief Adds a speciation rate for application post-simulation.
         * @param speciation_rate_in the speciation rate to add.
         */
        void addSpeciationRate(const double &speciation_rate_in);

        /**
         * @brief Adds metacommunity parameters for applying post-simulation.
         * @param metacommunity_size the number of individuals in the metacommunity
         * @param metacommunity_speciation_rate the speciation rate for the metacommunity
         * @param metacommunity_option the metacommunity option
         * @param metacommunity_reference the external database community reference
         */
        void addMetacommunityParameters(const unsigned long &metacommunity_size,
                                        const double &metacommunity_speciation_rate,
                                        const string &metacommunity_option,
                                        const unsigned long &metacommunity_reference);

        /**
         * @brief Adds the protracted parameters for applying post-simulation.
         * @param min_speciation_gen the minimum number of generations before speciation is permitted
         * @param max_speciation_gen the maximum number of generations a lineage can persist for
         */
        void addProtractedParameters(const double &min_speciation_gen, const double &max_speciation_gen);

        /**
         * @brief Sets the deme size and sampling proportion for the simulation.
         * @param deme_in the number of individuals existing in the landscape
         * @param deme_sample_in the proportion of individuals to sample
         */
        void setDeme(const unsigned long &deme_in, const double &deme_sample_in);

        /**
         * @brief Gets the sql database.
         * @return path to the sql database
         */
        string getSQLDatabase();

        /**
         * @brief Gets the output directory for the simulation.
         * @return the path to the output directory
         */
        string getOutputDirectory();

        /**
         * @brief Getter for the simulation seed.
         * @return Returns the seed
         */
        long long getSeed() override;

        /**
         * @brief Gets the job type for the simulation.
         * This is a reference number for the jobs.
         * @return Returns the job type
         */
        long long getJobType() override;

        /**
         * @brief CLoses the connection to the SQL file.
         *
         * Together with resumeSQLConnection(), the functions wrap  all SQL calls from R, opening or closing the database.
         */
        void pauseSQLConnection();

        /**
         * @brief Opens the connection to the SQL file.
         *
         * Together with pauseSQLConnection(), the functions wrap  all SQL calls from R, opening or closing the database.
         */
        void resumeSQLConnection();

        /**
         * @brief Applies the provided speciation parameters to the coalescence tree.
         * Records data for all times stored in the reference_times vector.
         * @param file_in the file path to save the database to
         * @param use_spatial_in if true, records all spatial locations of lineages
         * @param sample_file a map containing spatial sampling data
         * @param use_fragments_in if true, detects fragments from the map, if a file path, uses the file coordinates
         * @param speciation_rates a vector of speciation rates to apply
         * @param times_list a vector of times to apply at
         * @param min_speciation_gen_in the minimum number of generations required for speciation in protracted simulations
         * @param max_speciation_gen_in the maximum number of generations required before speciation in protracted sims
         * @param metacommunity_size_in the metacommunity size for protracted simulations
         * @param metacommunity_speciation_rate_in the metacommunity speciation rate for protracted simulations
         */
        void applySpeciation(const string &file_in,
                             const bool &use_spatial_in,
                             const string &sample_file,
                             const string &use_fragments_in,
                             vector<double> times_list);

        /**
         * @brief Calculates the abundance of each species and returns a dataframe containing species ids and abundances.
         * @return Dataframe containing species ids and abundances
         */
        Rcpp::DataFrame getSpeciesAbundances(const unsigned long &community_reference);

        /**
         * @brief Gets the number of species in the most recent calculation.
         * @return the number of species in the most recently-calculated coalescence tree
         */
        unsigned long getSpeciesRichness(const unsigned long &community_reference);

        /**
         * @brief Gets the lat calculated species richness from the community.
         * @return the number of species
         */
        unsigned long getLastSpeciesRichness();

        /**
         * @brief Checks that the main simulation has written to the in-memory database.
         */
        void checkWrittenMainSim();

        /**
         * @brief Writes the output to the previously specified database.
         */
        void output();

        /**
         * @brief Sets the logging mode to true or false.
         * @param log_mode if true, logs all messages to console
         */
        void setLoggingMode(bool log_mode);

        /**
         * @brief Sets the multiple_output variable - note that this doesn't override the multiple output variable if it
         * has already been set.
         * @param multiple_output
         */
        void setMultipleOutput(const bool &multiple_output);

        /**
         * @brief Gets the value of the multiple_output variable
         * @return true or false whether the simulation has applied multiple speciation rates
         * (and therefore needs to output)
         */
        bool getMultipleOutput();

        /**
         * @brief Applies the provided speciation parameters to the completed simulation.
         * @param specSimParameters
         */
        void apply(shared_ptr<SpecSimParameters> specSimParameters);

        /**
         * @brief Calls Tree::setup() to act as a wrapper accessible by R without extra classes.
         */
        void setupR() override;

        /**
         * @brief Calls Tree::runSimulation() to act as a wrapper accessible by R without extra classes.
         * @return bool true if simulation completes successfully
         */
        bool runSimulation();

    };
}
#endif //RCOALESCENCE_RTREE_H
