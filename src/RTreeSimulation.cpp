// This file is part of rcoalescence project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.
/**
 * @author Sam Thompson
 * @file RTree.cpp
 * @brief R wrapper for running non-spatial neutral models.
 *
 * Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */
#include "RTreeSimulation.h"
#include <utility>

RTreeSimulation::RTreeSimulation()
		: Tree(), spec_sim_parameters(make_shared<SpecSimParameters>()), multiple_output(false)
{
	output_database = "not_set";
	setLoggingMode(false);
}

RTreeSimulation::~RTreeSimulation()
{
	community.destroyDbConnection();
}

void RTreeSimulation::setKeyParameters(const long long &task_in, const long long &seed_in,
									   const string &output_directory_in,
									   const unsigned long &max_time_in, const unsigned long &desired_specnum_in,
									   vector<double> times_list)
{
	sim_parameters->setKeyParameters(task_in, seed_in, output_directory_in, max_time_in, desired_specnum_in,
									 "null");
	if(!times_list.empty())
	{
		for(auto &times: times_list)
		{
			sim_parameters->times.push_back(times);
		}
		times_file = "set";
		sim_parameters->times_file = "set";
	}
}

void RTreeSimulation::setMinSpeciationRate(const long double &spec_in)
{
	sim_parameters->setSpeciationParameters(spec_in, sim_parameters->is_protracted, sim_parameters->min_speciation_gen,
											sim_parameters->max_speciation_gen);
}

void RTreeSimulation::setSimulationProtractedParameters(const bool &protracted_in, const double &min_speciation_gen_in,
														const double &max_speciation_gen_in)
{
	sim_parameters->setSpeciationParameters(sim_parameters->spec,
											protracted_in, min_speciation_gen_in, max_speciation_gen_in);
}

void RTreeSimulation::addSpeciationRate(const double &speciation_rate_in)
{
	spec_sim_parameters->all_speciation_rates.insert(speciation_rate_in);
}

void RTreeSimulation::addMetacommunityParameters(const unsigned long &metacommunity_size,
												 const double &metacommunity_speciation_rate,
												 const string &metacommunity_option,
												 const unsigned long &metacommunity_reference)
{
	spec_sim_parameters->addMetacommunityParameters(metacommunity_size, metacommunity_speciation_rate,
													metacommunity_option, metacommunity_reference);
}

void RTreeSimulation::addProtractedParameters(const double &min_speciation_gen, const double &max_speciation_gen)
{
	spec_sim_parameters->addProtractedParameters(min_speciation_gen, max_speciation_gen);
}

void RTreeSimulation::setDeme(const unsigned long &deme_in, const double &deme_sample_in)
{
	sim_parameters->deme = deme_in;
	sim_parameters->deme_sample = deme_sample_in;
}

string RTreeSimulation::getSQLDatabase()
{
	return sql_output_database;
}

void RTreeSimulation::apply(shared_ptr<SpecSimParameters> specSimParameters)
{
	if(!sim_complete)
	{
		throw FatalException("Cannot apply speciation rates to an incomplete simulation.");
	}
	else
	{
		community.setSimParameters(sim_parameters);
		community.doApplicationInternal(std::move(specSimParameters), data);
	}
}

void RTreeSimulation::applySpeciation(const string &file_in, const bool &use_spatial_in, const string &sample_file,
									  const string &use_fragments_in, vector<double> times_list)
{
	checkDatabaseSet();
	spec_sim_parameters->setup(file_in, use_spatial_in, sample_file, times_list, use_fragments_in);
	apply(spec_sim_parameters);
	spec_sim_parameters->wipe();
}

Rcpp::DataFrame RTreeSimulation::getSpeciesAbundances(const unsigned long &community_reference)
{
	auto row = community.getSpeciesAbundances(community_reference);
	Rcpp::IntegerVector species_ids(row->size());
	Rcpp::IntegerVector no_individuals(row->size());
	unsigned long i = 0;
	for(auto const &x : *row)
	{
		species_ids[i] = x.first;
		no_individuals[i] = x.second;
	}
	Rcpp::DataFrame out_df = Rcpp::DataFrame::create(Rcpp::Named("species_id") = species_ids,
													 Rcpp::Named("no_individuals") = no_individuals);
	return out_df;
}

unsigned long RTreeSimulation::getSpeciesRichness(const unsigned long &community_reference)
{
	return community.getSpeciesRichness(community_reference);
}

unsigned long RTreeSimulation::getLastSpeciesRichness()
{
	return community.getSpeciesNumber();
}

void RTreeSimulation::checkDatabaseSet()
{
	if(!community.isSetDatabase())
	{
		openSQLDatabase();
		community.setDatabase(database);
	}
}

void RTreeSimulation::output()
{
	setupOutputDirectory();
	checkDatabaseSet();
	spec_sim_parameters->filename = sql_output_database;
	community.output();
	community.destroyDbConnection();
}

void RTreeSimulation::setLoggingMode(bool log_mode)
{
	logging_mode = log_mode;
}

void RTreeSimulation::setMultipleOutput(const bool &multiple_output)
{
	if(!this->multiple_output)
	{
		this->multiple_output = multiple_output;
	}
}

bool RTreeSimulation::getMultipleOutput()
{
	return multiple_output;
}

void RTreeSimulation::setup()
{
	Tree::setup();
}

bool RTreeSimulation::runSimulation()
{
	bool completed = Tree::runSimulation();
}

