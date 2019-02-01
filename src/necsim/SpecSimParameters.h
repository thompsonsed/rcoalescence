// This file is part of necsim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @author Sam Thompson
 * @file SpecSimParameters.h
 * @brief Contains parameters for applying speciation rates post-simulation.
 *
 * Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */



#ifndef SPECIATIONCOUNTER_SPECSIMPARAMETERS_H
#define SPECIATIONCOUNTER_SPECSIMPARAMETERS_H

#include <string>
#include <utility>
#include <vector>
#include "ConfigParser.h"
#include "custom_exceptions.h"
#include "double_comparison.h"
#include "parameters.h"

using namespace std;

/**
 * @class SpecSimParameters
 * @brief Contains the simulation parameters that are read from the command line.
 *
 */
struct SpecSimParameters
{
    bool use_spatial;
    bool bMultiRun;
    bool use_fragments;
    string filename;
    set<double> all_speciation_rates;
    string samplemask;
    string times_file;
    set<double> all_times;
    string fragment_config_file;
    vector<ProtractedSpeciationParameters> protracted_parameters;
    MetacommunitiesArray metacommunity_parameters;

    SpecSimParameters() : use_spatial(false), bMultiRun(false), use_fragments(false), filename("none"),
                          all_speciation_rates(), samplemask("none"), times_file("null"), all_times(),
                          fragment_config_file("none"), protracted_parameters(), metacommunity_parameters()
    {

    }

    SpecSimParameters(const string &fragment_config_file) : fragment_config_file(fragment_config_file){}

    /**
     * @brief Adds a speciation rate to the speciation parameters list.
     * @param speciation_rate the speciation rate to add
     */
    void addSpeciationRate(const double &speciation_rate)
    {
        all_speciation_rates.insert(speciation_rate);
        bMultiRun = all_speciation_rates.size() > 1;
    }

    /**
     * @brief Sets the application arguments for the inputs.
     * @param file_in the database to apply speciation rates to
     * @param use_spatial_in if true, record full spatial data
     * @param sample_file the sample file to select lineages from the map
     * @param use_fragments_in fragment file, or "T"/"F" for automatic detection/no detection
     */
    void setup(string file_in, bool use_spatial_in, string sample_file, const string &use_fragments_in)
    {
        filename = std::move(file_in);
        use_spatial = use_spatial_in;
        samplemask = std::move(sample_file);
        use_fragments = !(use_fragments_in == "F");
        fragment_config_file = use_fragments_in;
    }

    /**
     * @brief Sets the application arguments for the inputs.
     * @param file_in the database to apply speciation rates to
     * @param use_spatial_in if true, record full spatial data
     * @param sample_file the sample file to select lineages from the map
     * @param times vector of times to apply
     * @param use_fragments_in fragment file, or "T"/"F" for automatic detection/no detection
     * @param speciation_rates the speciation rates to apply
     */
    void setup(string file_in, bool use_spatial_in, string sample_file, const vector<double> &times,
               const string &use_fragments_in, vector<double> speciation_rates)
    {
        setup(file_in, use_spatial_in, sample_file, times, use_fragments_in);
        for(const auto &speciation_rate : speciation_rates)
        {
            addSpeciationRate(speciation_rate);
        }
    }

    /**
     * @brief Sets the application arguments for the inputs. Overloaded version without speciation rates.
     * @param file_in the database to apply speciation rates to
     * @param use_spatial_in if true, record full spatial data
     * @param sample_file the sample file to select lineages from the map
     * @param times vector of times to apply
     * @param use_fragments_in fragment file, or "T"/"F" for automatic detection/no detection
     */
    void setup(string file_in, bool use_spatial_in, string sample_file, const vector<double> &times,
               const string &use_fragments_in)
    {
        setup(file_in, use_spatial_in, sample_file, use_fragments_in);
        if(times.empty() && all_times.empty())
        {
            times_file = "null";
            all_times.insert(0.0);
        }
        else
        {
            for(const auto &time: times)
            {
                addTime(time);
            }
        }
    }



    /**
     * @brief Sets the metacommunity parameters for the simulation.
     * @param metacommunity_size_in the number of individuals in the metacommunity
     * @param metacommunity_speciation_rate_in the speciation rate for the metacommunity
     */
    void addMetacommunityParameters(const unsigned long &metacommunity_size_in,
                                    const double &metacommunity_speciation_rate_in,
                                    const string &metacommunity_option_in,
                                    const unsigned long &metacommunity_reference_in)
    {
        MetacommunityParameters tmp_meta_parameters = MetacommunityParameters(0, metacommunity_size_in,
                                                                              metacommunity_speciation_rate_in,
                                                                              metacommunity_option_in,
                                                                              metacommunity_reference_in);
        if(!metacommunity_parameters.hasOption(tmp_meta_parameters))
        {
            metacommunity_parameters.addNew(tmp_meta_parameters);
        }
        else
        {
            stringstream ss;
            ss << "Parameters already added for metacommunity with size " << metacommunity_size_in << ", ";
            ss << "speciation rate " << metacommunity_speciation_rate_in << ", " << "using option ";
            ss << metacommunity_option_in << " and reference " << metacommunity_reference_in << endl;
            writeInfo(ss.str());
        }
    }

    /**
     * @brief Import the time config file, if there is one
     */
    void importTimeConfig()
    {
        if(times_file == "null")
        {
            all_times.insert(0.0);
        }
        else
        {
            vector<string> tmpimport;
            ConfigParser tmpconfig;
            tmpconfig.setConfig(times_file, false);
            tmpconfig.importConfig(tmpimport);
            for(const auto &i : tmpimport)
            {
                all_times.insert(stod(i));
            }
        }
    }

    /**
     * @brief Deletes all the parameters.
     */
    void wipe()
    {
        use_spatial = false;
        bMultiRun = false;
        use_fragments = false;
        filename = "";
        all_speciation_rates.clear();
        samplemask = "";
        times_file = "";
        all_times.clear();
        fragment_config_file = "";
        protracted_parameters.clear();
        metacommunity_parameters.clear();
    }

    /**
     * @brief Adds an additional time to the times vector.
     * @param time a time to calculate speciation rates at
     */
    void addTime(double time)
    {
        all_times.insert(time);
        if(all_times.size() > 1)
        {
            times_file = "set";
        }
    }

    /**
     * @brief Adds a set of protracted speciation parameters to the protracted parameters vector
     * @param proc_spec_min the minimum protracted speciation generation
     * @param proc_spec_max the maximum protracted speciation generation
     */
    void addProtractedParameters(double proc_spec_min, double proc_spec_max)
    {
        ProtractedSpeciationParameters tmp;
        tmp.min_speciation_gen = proc_spec_min;
        tmp.max_speciation_gen = proc_spec_max;
        protracted_parameters.emplace_back(tmp);
    }
};

#endif //SPECIATIONCOUNTER_SPECSIMPARAMETERS_H
