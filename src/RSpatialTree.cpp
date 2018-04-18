//This file is part of rcoalescence project which is released under BSD-3 license.
//Visit https://opensource.org/licenses/BSD-3-Clause) for full license details.
#include "RSpatialTree.h"
#include <utility> 
#include "RLogging.h"

bool logging_mode = false;

RSpatialTree::RSpatialTree() : SpatialTree()
{
	output_database = "not_set";
	setLoggingMode(false);
}


RSpatialTree::~RSpatialTree()
= default;

void RSpatialTree::setKeyParameters(const long long &task_in, const long long &seed_in,
									const string &output_directory_in,
									const unsigned long &max_time_in, const unsigned long &desired_specnum_in,
									vector<double> times_list)
{
	sim_parameters.setKeyParameters(task_in, seed_in, output_directory_in, max_time_in, desired_specnum_in,
									"null");
	if(!times_list.empty())
	{
		for(auto &times: times_list)
		{
			sim_parameters.times.push_back(times);
		}
		times_file = "set";
		sim_parameters.times_file = "set";
	}
}

void RSpatialTree::setSpeciationParameters(const long double &spec_in, bool is_protracted_in,
											 const double &min_speciation_gen_in, const double &max_speciation_gen_in)
{
	sim_parameters.setSpeciationParameters(spec_in, is_protracted_in, min_speciation_gen_in, max_speciation_gen_in);
}

void RSpatialTree::setDispersalParameters(const double &sigma_in, const string &dispersal_method_in, const double &tau_in,
											const double &m_prob_in, const double &cutoff_in,
											const double &dispersal_relative_cost_in, bool restrict_self_in,
											const string &landscape_type_in, const string &dispersal_file_in,
											const string &reproduction_file_in)
{
	sim_parameters.setDispersalParameters(dispersal_method_in, sigma_in, tau_in, m_prob_in, cutoff_in,
										  dispersal_relative_cost_in, restrict_self_in, landscape_type_in,
										  dispersal_file_in, reproduction_file_in);
}

void RSpatialTree::setHistoricalMapParameters(const string &historical_fine_file_map_in,
											  const string &historical_coarse_map_file_in,
											  const double &gen_since_historical_in,
											  const double &habitat_change_rate_in)
{
	sim_parameters.setHistoricalMapParameters(historical_fine_file_map_in, historical_coarse_map_file_in,
											gen_since_historical_in, habitat_change_rate_in);
}

void RSpatialTree::addHistoricalMap(const string &fine_map, const string &coarse_map, const double &time,
									const double &rate)
{
	if(sim_parameters.historical_fine_map_file == "none")
	{
		sim_parameters.historical_fine_map_file = fine_map;
		sim_parameters.historical_coarse_map_file = coarse_map;
	}
	paths_fine.push_back(fine_map);
	paths_coarse.push_back(coarse_map);
	unsigned long last_number = 0;
	if(!numbers_fine.empty())
	{
		last_number = numbers_fine.back() + 1;
	}
	numbers_fine.push_back(last_number);
	numbers_coarse.push_back(last_number);
	times_fine.push_back(time);
	times_coarse.push_back(time);
	rates_fine.push_back(rate);
	rates_coarse.push_back(rate);
}

void RSpatialTree::setMapParameters(const string &fine_map_file_in, const string &coarse_map_file_in,
									  const string &sample_mask_file_in, const unsigned long &grid_x_size_in,
									  const unsigned long &grid_y_size_in, const unsigned long &sample_x_size_in,
									  const unsigned long &sample_y_size_in, const unsigned long &sample_x_offset_in,
									  const unsigned long &sample_y_offset_in, const unsigned long &fine_map_x_size_in,
									  const unsigned long &fine_map_y_size_in,
									  const unsigned long &fine_map_x_offset_in,
									  const unsigned long &fine_map_y_offset_in,
									  const unsigned long &coarse_map_x_size_in,
									  const unsigned long &coarse_map_y_size_in,
									  const unsigned long &coarse_map_x_offset_in,
									  const unsigned long &coarse_map_y_offset_in,
									  const unsigned long &coarse_map_scale_in, const unsigned long &deme_in,
									  const double &deme_sample_in, bool uses_spatial_sampling_in)
{
	sim_parameters.setMapParameters(fine_map_file_in, coarse_map_file_in, sample_mask_file_in, grid_x_size_in,
									grid_y_size_in, sample_x_size_in, sample_y_size_in, sample_x_offset_in,
									sample_y_offset_in, fine_map_x_size_in, fine_map_y_size_in, fine_map_x_offset_in,
									fine_map_y_offset_in, coarse_map_x_size_in, coarse_map_y_size_in,
									coarse_map_x_offset_in, coarse_map_y_offset_in, coarse_map_scale_in,
									deme_in, deme_sample_in, uses_spatial_sampling_in);
}

string RSpatialTree::getSQLDatabase()
{
	return sql_output_database;
}

void RSpatialTree::setup()
{
	if(!paths_fine.empty())
	{
		sim_parameters.setHistoricalMapParameters(paths_fine, numbers_fine, rates_fine, times_fine, paths_coarse,
												  numbers_coarse, rates_coarse, times_coarse);
	}
	SpatialTree::setup();
}

bool RSpatialTree::runSimulation()
{
	if(!paths_fine.empty())
	{
		sim_parameters.setHistoricalMapParameters(paths_fine, numbers_fine, rates_fine, times_fine, paths_coarse,
												  numbers_coarse, rates_coarse, times_coarse);
		sim_parameters.setHistorical(0);
		landscape.resetHistorical();
	}
	return Tree::runSimulation();
}

void RSpatialTree::apply(SpecSimParameters *specSimParameters)
{
	if(!sim_complete)
	{
		throw FatalException("Cannot apply speciation rates to an incomplete simulation.");
	}
	else
	{
		community.setSimParameters(&sim_parameters);
		community.doApplicationInternal(specSimParameters, &data);
	}
}

void RSpatialTree::applySpeciation(const string &file_in, const bool &use_spatial_in, const string &sample_file,
								   const string &use_fragments_in, vector<double> speciation_rates,
								   vector<double> times_list,
								   const double &min_speciation_gen_in, const double &max_speciation_gen_in,
								   const unsigned long &metacommunity_size_in,
								   const double &metacommunity_speciation_rate_in)
{
	checkDatabaseSet();
	spec_sim_parameters.wipe();
	spec_sim_parameters.setup(file_in, use_spatial_in, sample_file, std::move(times_list),
							  use_fragments_in, std::move(speciation_rates),
							  min_speciation_gen_in, max_speciation_gen_in,
							  metacommunity_size_in, metacommunity_speciation_rate_in);
	apply(&spec_sim_parameters);
}

Rcpp::DataFrame RSpatialTree::getSpeciesAbundances()
{
	auto row = community.getRowOut();
	Rcpp::IntegerVector species_ids(row.size());
	Rcpp::IntegerVector no_individuals(row.size());
	for(unsigned long i = 0; i < row.size(); i ++)
	{
		species_ids[i] = i;
		no_individuals[i] = row[i];
	}
	Rcpp::DataFrame out_df = Rcpp::DataFrame::create(Rcpp::Named("species_id")=species_ids,
													 Rcpp::Named("no_individuals")=no_individuals);
	return out_df;
}

double RSpatialTree::getSpeciesRichness()
{
	return community.getSpeciesNumber();
}

void RSpatialTree::checkDatabaseSet()
{
	if(!community.isSetDatabase())
	{
		openSQLDatabase();
		community.setDatabase(database);
	}
}

void RSpatialTree::output()
{
	sortData();
	sqlCreate();
//	outputData();
	checkDatabaseSet();
	spec_sim_parameters.filename = sql_output_database;
	community.output();
}

void RSpatialTree::setLoggingMode(bool log_mode)
{
	logging_mode = log_mode;
}

RCPP_EXPOSED_CLASS(RSpatialTree);
RCPP_MODULE(coalescenceModule) {
	using namespace Rcpp;
	class_<RSpatialTree>("RSpatialTree", "Simulates spatially-explicit neutral models.")
			.constructor("initialises the spatial tree")
			.method("._setKeyParameters", &RSpatialTree::setKeyParameters,
					"Sets the main simulation parameters for house-keeping purposes.")
			.method("._setSpeciationParameters", &RSpatialTree::setSpeciationParameters,
					"Sets the speciation parameters for the simulation.")
			.method("._setDispersalParameters", &RSpatialTree::setDispersalParameters,
					"Sets the dispersal parameters for the simulation.")
			.method("._setHistoricalMapParameters", &RSpatialTree::setHistoricalMapParameters,
					"Sets the historical map parameters for the simulation.")
			.method("._setMapParameters", &RSpatialTree::setMapParameters,
					"Sets the map parameters for the simulation, including dimensions and offsets.")
			.method("setup", &RSpatialTree::setup, "Completes all set-up routines for the simulation, "
					"including assigning memory and importing the required files.")
			.method("._addHistoricalMap", &RSpatialTree::addHistoricalMap, "Adds a historical map with the desired parameter "
					"set")
			.method("runSimulation", &RSpatialTree::runSimulation, "Performs the actual simulation, "
					"returning true if the simulation is complete.")
			.method("._applySpeciationRates", &RSpatialTree::applySpeciation, "Applies the provided speciation parameters"
					"to the completed simulation.")
			.method("getSpeciesAbundances", &RSpatialTree::getSpeciesAbundances, "Gets the species abundances for the "
					"last calculated set of speciation parameters.")
			.method("._getSpeciesRichness", &RSpatialTree::getSpeciesRichness, "Gets the species richness for the last "
					"calculated set of speciation parameters.")
			.method("._output", &RSpatialTree::output, "Writes the simulation to the output database.")
			.method("setLoggingMode", &RSpatialTree::setLoggingMode, "Sets the logging mode for the simulation object.")
			.method("._getSQLDatabase", &RSpatialTree::getSQLDatabase)
			.field("output_database", &RSpatialTree::output_database)
			;
}

