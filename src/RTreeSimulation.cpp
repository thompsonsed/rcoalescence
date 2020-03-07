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
using namespace rcoalescence;
using namespace std;
using namespace necsim;
RTreeSimulation::RTreeSimulation()
        : Tree(), spec_sim_parameters(make_shared<SpecSimParameters>()), metacommunity(), multiple_output(false),
          has_outputted(false), has_written_main_sim(false), uses_metacommunity(false)
{
    output_database = "not_set";
    setLoggingMode(false);
    community.setInternalDatabase();
    metacommunity.setInternalDatabase();
}

RTreeSimulation::~RTreeSimulation()
{
    community.closeSqlConnection();
    metacommunity.closeSqlConnection();
}

void RTreeSimulation::setKeyParameters(const long long &job_type, const long long &seed_in,
                                       const string &output_directory_in,
                                       const unsigned long &max_time_in, const unsigned long &desired_specnum_in,
                                       vector<double> times_list)
{
    sim_parameters->setKeyParameters(job_type, seed_in, output_directory_in, max_time_in, desired_specnum_in,
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
    stringstream ss;
    ss << "Adding metacommunity parameters with size = " << metacommunity_size << ", speciation rate = ";
    ss << metacommunity_speciation_rate << ", option = " << metacommunity_option << " and external reference = ";
    ss << metacommunity_reference << endl;
    writeInfo(ss.str());
    spec_sim_parameters->addMetacommunityParameters(metacommunity_size, metacommunity_speciation_rate,
                                                    metacommunity_option, metacommunity_reference);
}

void RTreeSimulation::addProtractedParameters(const double &min_speciation_gen, const double &max_speciation_gen)
{
    stringstream ss;
    ss << "Adding protracted parameters with min = " << min_speciation_gen << ", max = " << max_speciation_gen << endl;
    writeInfo(ss.str());
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

string RTreeSimulation::getOutputDirectory()
{
    return out_directory;
}

long long RTreeSimulation::getSeed()
{
    return seed;
}

long long RTreeSimulation::getJobType()
{
    return job_type;
}

void RTreeSimulation::apply(shared_ptr<SpecSimParameters> specSimParameters)
{
    if(!sim_complete)
    {
        throw FatalException("Cannot apply speciation rates to an incomplete simulation.");
    }
    else
    {
        if(specSimParameters->metacommunity_parameters.hasMetacommunityOption() || uses_metacommunity)
        {
            if(!uses_metacommunity)
            {
                metacommunity.setSpecSimParameters(spec_sim_parameters);
                metacommunity.setupInternal(sim_parameters, database);
            }
            metacommunity.applyNoOutput(std::move(specSimParameters), data);
            uses_metacommunity = true;
        }
        else
        {
            community.setSimParameters(sim_parameters);
            community.doApplicationInternal(std::move(specSimParameters), data);
        }
    }
}

void RTreeSimulation::applySpeciation(const string &file_in, const bool &use_spatial_in, const string &sample_file,
                                      const string &use_fragments_in, vector<double> times_list)
{
    checkWrittenMainSim();
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
    if(uses_metacommunity)
    {
        return metacommunity.getSpeciesNumber();
    }
    return community.getSpeciesNumber();
}


void RTreeSimulation::checkWrittenMainSim()
{
    if(!has_written_main_sim)
    {
        sortData();
        sqlCreate();
        community.setSpecSimParameters(spec_sim_parameters);
        setupCommunity();
        has_written_main_sim = true;
    }
}

void RTreeSimulation::output()
{
    if(has_outputted)
    {
        throw FatalException("Output database has already been generated.");
    }
    setupOutputDirectory();
    checkWrittenMainSim();
    spec_sim_parameters->filename = sql_output_database;
    // TODO move this to debug or remove
    if(database == nullptr || community.isDatabaseNullPtr())
    {
        throw FatalException("Database pointer is null before output. Please report this bug.");
    }
    if(uses_metacommunity)
    {
        metacommunity.output();
        metacommunity.closeSqlConnection();
    }
    else
    {
        community.output();
    }
    community.closeSqlConnection();
    output_database = sql_output_database;
    has_outputted = true;
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
    return Tree::runSimulation();
}

