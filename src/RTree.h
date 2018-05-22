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

#include <Rcpp.h>
//#include <Rcpp/r/headers.h>
//#include <Rcpp/DataFrame.h>
#include <necsim/SpecSimParameters.h>
#include <necsim/Tree.h>
#include "RLogging.h"

using namespace std;
class RTree : public virtual Tree
{
protected:

	SpecSimParameters spec_sim_parameters{};
public:
	string output_database;
	RTree();
	~RTree() override;


	/**
	 * @brief Sets the main simulation parameters in the sim_parameters object.
	 * @param task_in the task reference number, used for file referencing
	 * @param seed_in the seed to set random number generation
	 * @param output_directory_in the output directory
	 * @param max_time_in the maximum time to simulate for
	 * @param desired_specnum_in the desired number of species to aim towards (currently not functional)
	 * @param times_list the file containing a list of temporal sampling points
	 */
	void setKeyParameters(const long long &task_in, const long long &seed_in, const string &output_directory_in,
						  const unsigned long &max_time_in, const unsigned long &desired_specnum_in,
						  vector<double> times_list);

	/**
	 * @brief Sets the speciation parameters for the simulation in the sim_parameters object.
	 * @param spec_in the speciation rate to use
	 * @param is_protracted_in if true, simulates as a protracted simulation
	 * @param min_speciation_gen_in the minimum speciation generation for protracted simulations
	 * @param max_speciation_gen_in the maximum speciation generation for protracted simulations
	 */
	void setSpeciationParameters(const long double &spec_in, bool is_protracted_in, const double &min_speciation_gen_in,
								 const double &max_speciation_gen_in);

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
	void applySpeciation(const string &file_in, const bool &use_spatial_in, const string &sample_file,
						 const string &use_fragments_in, vector<double> speciation_rates, vector<double> times_list,
						 vector<double> min_speciation_gen_in, vector<double> max_speciation_gen_in,
						 const unsigned long &metacommunity_size_in,
						 const double &metacommunity_speciation_rate_in);

	/**
	 * @brief Calculates the abundance of each species and returns a dataframe containing species ids and abundances.
	 * @return Dataframe containing species ids and abundances
	 */
	Rcpp::DataFrame getSpeciesAbundances();


	/**
	 * @brief Gets the number of species in the most recent calculation.
	 * @return the number of species in the most recently-calculated coalescence tree
	 */
	unsigned long getSpeciesRichness();

	/**
	 * @brief Ensures that a connection is made to the output database.
	 */
	void checkDatabaseSet();

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
	 * @brief Applies the provided speciation parameters to the completed simulation.
	 * @param specSimParameters
	 */
	void apply(SpecSimParameters *specSimParameters);

	/**
	 * @brief Calls Tree::setup() to act as a wrapper accessible by R without extra classes.
	 */
	void setup() override;

	/**
	 * @brief Calls Tree::runSimulation() to act as a wrapper accessible by R without extra classes.
	 * @return bool true if simulation completes successfully
	 */
	bool runSimulation() override;

};

#endif //RCOALESCENCE_RTREE_H
