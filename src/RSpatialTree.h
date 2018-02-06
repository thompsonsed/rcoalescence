//
// Created by Sam Thompson on 04/02/2018.
//

#ifndef RCOALESCENCE_R_SPATIAL_TREE_H
#define RCOALESCENCE_R_SPATIAL_TREE_H
#include "Rcpp.h"
#include <Rcpp/DataFrame.h>
#include <Rcpp/vector/instantiation.h>
#include "necsim/SpatialTree.h"
//' @export
//' Class containing routines for spatially-explicit neutral models performed on the provided maps.
class RSpatialTree : public SpatialTree
{
public:
	RSpatialTree();
	~RSpatialTree() override;

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
	 * @brief Sets the dispersal parameters for the simulation.
 	 * @param sigma_in the sigma value for a normal distribution
	 * @param dispersal_method_in the method of individuals dispersing (normal, fat-tailed or norm-uniform)
	 * @param tau_in the tau value for the fat-tailed distribution
	 * @param m_prob_in the probability of uniform dispersal for the norm-uniform distribution
	 * @param cutoff_in the maximum dispersal distance for the uniform distribution
	 * @param dispersal_relative_cost_in the relative cost of dispersing through non-forest
	 * @param restrict_self_in if true, prevents dispersal from the same cell
	 * @param landscape_type_in the landscape type (infinite, tiled or closed)
	 * @param dispersal_file_in a map of dispersal probabilities
	 * @param reproduction_file_in a map of reproduction probabilities
	 */
	void setDispersalParameters( const double &sigma_in, const string &dispersal_method_in, const double &tau_in,
								const double &m_prob_in, const double &cutoff_in,
								const double &dispersal_relative_cost_in, bool restrict_self_in,
								const string &landscape_type_in, const string &dispersal_file_in,
								const string &reproduction_file_in);

	/**
	 * @brief Sets the pristine map parameters for the simulation.
	 * @param pristine_fine_file_map_in the fine resolution pristine file
	 * @param pristine_coarse_map_file_in the coarse resolution pristine file
	 * @param gen_since_pristine_in the number of generations since the pristine state was achieved
	 * @param habitat_change_rate_in the rate of habitat change towards the pristine state
	 */
	void setPristineMapParameters(const string &pristine_fine_file_map_in, const string &pristine_coarse_map_file_in,
								  const double &gen_since_pristine_in, const double &habitat_change_rate_in);

	/**
	 * @brief Sets the map parameters for the simulation.
	 * @param fine_map_file_in the fine resolution density map
	 * @param coarse_map_file_in the coarse resolution density map
	 * @param sample_mask_file_in the spatial sampling mask
	 * @param grid_x_size_in the x dimension of the grid
	 * @param grid_y_size_in the y dimension of the grid
	 * @param sample_x_size_in the x dimension of the sample mask
	 * @param sample_y_size_in the y dimension of the sample mask
	 * @param sample_x_offset_in the x offset of the sample mask from the grid
	 * @param sample_y_offset_in the y offset of the sample mask from the grid
	 * @param fine_map_x_size_in the x dimension of the fine map
	 * @param fine_map_y_size_in the y dimension of the fine map
	 * @param fine_map_x_offset_in the x offset of the fine map from the sample mask
	 * @param fine_map_y_offset_in the y offset of the fine map from the sample mask
	 * @param coarse_map_x_size_in the x dimension of the coarse map
	 * @param coarse_map_y_size_in the y dimension of the coarse map
	 * @param coarse_map_x_offset_in the x offset of the coarse map from the fine map
	 * @param coarse_map_y_offset_in the y offset of the coarse map from the fine map
	 * @param coarse_map_scale_in the scale of the coarse map compared to the fine map
	 * @param deme_in the number of individuals per cell
	 * @param deme_sample_in the proportion of individuals to sample from each cell
	 * @param uses_spatial_sampling_in if the sample mask denotes differing spatial sampling proportions
	 */
	void setMapParameters(const string &fine_map_file_in, const string &coarse_map_file_in,
						  const string &sample_mask_file_in, const unsigned long &grid_x_size_in,
						  const unsigned long &grid_y_size_in, const unsigned long &sample_x_size_in,
						  const unsigned long &sample_y_size_in, const unsigned long &sample_x_offset_in,
						  const unsigned long &sample_y_offset_in, const unsigned long &fine_map_x_size_in,
						  const unsigned long &fine_map_y_size_in, const unsigned long &fine_map_x_offset_in,
						  const unsigned long &fine_map_y_offset_in, const unsigned long &coarse_map_x_size_in,
						  const unsigned long &coarse_map_y_size_in, const unsigned long &coarse_map_x_offset_in,
						  const unsigned long &coarse_map_y_offset_in, const unsigned long &coarse_map_scale_in,
						  const unsigned long &deme_in, const double &deme_sample_in, bool uses_spatial_sampling_in);

	/**
	 * @brief Calls SpatialTree::setup() to act as a wrapper accessible by R without extra classes.
	 */
	void setup() override;

	/**
	 * @brief Calls SpatialTree::runSimulation() to act as a wrapper accessible by R without extra classes.
	 * @return bool true if simulation completes successfully
	 */
	bool runSimulation() override;

	/**
	 * @brief Applies the provided speciation parameters to the completed simulation.
	 * @param specSimParameters
	 */
	void apply(SpecSimParameters *specSimParameters);

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
						 const double &min_speciation_gen_in, const double &max_speciation_gen_in,
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
	double getSpeciesRichness();

	/**
	 * @brief Writes the output to the previously specified database.
	 */
	void output();


};
#endif //RCOALESCENCE_R_SPATIAL_TREE_H
