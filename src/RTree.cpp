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
#include "RTree.h"

RTree::RTree() : Tree()
{
	output_database = "not_set";
	setLoggingMode(false);
}

RTree::~RTree()
= default;


void RTree::setKeyParameters(const long long &task_in, const long long &seed_in,
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

void RTree::setSpeciationParameters(const long double &spec_in, bool is_protracted_in,
										   const double &min_speciation_gen_in, const double &max_speciation_gen_in)
{
	sim_parameters.setSpeciationParameters(spec_in, is_protracted_in, min_speciation_gen_in, max_speciation_gen_in);
}

void RTree::setDeme(const unsigned long &deme_in, const double &deme_sample_in)
{
	sim_parameters.deme = deme_in;
	sim_parameters.deme_sample= deme_sample_in;
}

string RTree::getSQLDatabase()
{
	return sql_output_database;
}

void RTree::apply(SpecSimParameters *specSimParameters)
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

void RTree::applySpeciation(const string &file_in, const bool &use_spatial_in, const string &sample_file,
							const string &use_fragments_in, vector<double> speciation_rates,
							vector<double> times_list,
							vector<double> min_speciation_gen_in, vector<double> max_speciation_gen_in,
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

Rcpp::DataFrame RTree::getSpeciesAbundances()
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

unsigned long RTree::getSpeciesRichness()
{
	unsigned long species_number = 0;
	auto rows = community.getRowOut();
	for(unsigned long i = 0; i < rows.size(); i ++)
	{
		if(rows[i] > 0)
		{
			species_number ++;
		}
	}
	return species_number;
}

void RTree::checkDatabaseSet()
{
	if(!community.isSetDatabase())
	{
		openSQLDatabase();
		community.setDatabase(database);
	}
}

void RTree::output()
{
	sortData();
	sqlCreate();
//	outputData();
	checkDatabaseSet();
	spec_sim_parameters.filename = sql_output_database;
	community.output();
}


void RTree::setLoggingMode(bool log_mode)
{
}

void RTree::setup()
{
	Tree::setup();
}

bool RTree::runSimulation()
{
	return Tree::runSimulation();
}

