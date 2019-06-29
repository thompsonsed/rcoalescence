// This file is part of necsim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.
/**
 * @author Sam Thompson
 * @file SimParameters.h
 * @brief Stores and parses simulation parameters from the command line or a config file.
 *
 * Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */

#ifndef SPECIATIONCOUNTER_SIMPARAMETERS_H
#define SPECIATIONCOUNTER_SIMPARAMETERS_H

#include <string>
#include <utility>
#include <vector>
#include "ConfigParser.h"
#include "Logging.h"
#include "custom_exceptions.h"
#include "file_system.h"

using namespace std;
/************************************************************
					MAPVARS STRUCTURE
 ************************************************************/
/**
 * @struct SimParameters
 * @brief Stores and imports the variables required by the Map object.
 * Used to setting the Map variables in a more elegant way.
 */
struct SimParameters
{
    string fine_map_file, coarse_map_file, output_directory;
    string historical_fine_map_file, historical_coarse_map_file, sample_mask_file;
    // for file naming purposes.
    long long job_type{}, seed{};
    // the variables for the grid containing the initial individuals.
    unsigned long grid_x_size{}, grid_y_size{};
    // The variables for the sample grid, which may or may not be the same as the main simulation grid
    unsigned long sample_x_size{}, sample_y_size{};
    unsigned long sample_x_offset{}, sample_y_offset{};
    // The fine map variables at the same resolution as the grid.
    unsigned long fine_map_x_size{}, fine_map_y_size{}, fine_map_x_offset{}, fine_map_y_offset{};
    // the coarse map variables at a scaled resolution of the fine map.
    unsigned long coarse_map_x_size{}, coarse_map_y_size{}, coarse_map_x_offset{}, coarse_map_y_offset{};
    unsigned long coarse_map_scale{};
    unsigned long desired_specnum{};
    // the relative cost of moving through non-forest
    double dispersal_relative_cost{};
    // the size of each square of habitat in numbers of individuals
    double deme{};
    // the sample proportion,
    double deme_sample{};
    // the speciation rate.
    long double spec{0.0};
    // the variance of the dispersal kernel.
    double sigma{};
    // max time to run for
    unsigned long max_time{};
    // the number of generations since a historical landscape was encountered.
    double gen_since_historical{};
    // the transform rate of the forest from historical to modern forest.
    double habitat_change_rate{};
    // the fatness of the dispersal kernel
    double tau{};
    // dispersal method - should be one of [normal, fat-tail, norm-uniform]
    string dispersal_method;
    // the probability of selecting from a uniform dispersal kernel (for uniformally-modified dispersals)
    double m_prob{};
    // the cutoff for the normal dispersal in cells.
    double cutoff{};
    // if true, restricts dispersal from the same cell.
    bool restrict_self{};
    // file containing the points to record data from
    string times_file;
    // vector of times
    vector<double> times;
    // Stores the full list of configs imported from file
    ConfigParser configs;
    // Set to true if the oldest historical state has been reached.
    bool is_historical{};
    // if the sample file is not null, this variable tells us whether different points in space require different
    // numbers of individuals to be sampled. If this is the case, the actual values are read from the sample mask as a
    // proportion of individuals sampled, from 0-1. Otherwise, it is treated as a boolean mask, with values > 0.5
    // representing sampling in the cell.
    bool uses_spatial_sampling{};
    // This can be closed, infinite and tiled (which is also infinite)
    string landscape_type;
    // The protracted speciation current_metacommunity_parameters - these DON'T need to be stored upon pausing simulations
    bool is_protracted{};
    double min_speciation_gen{}, max_speciation_gen{};

    // a map of dispersal values, where each row corresponds to the probability of moving from one cell
    // to any other.
    string dispersal_file;

    // a map of relative death probabilities.
    string death_file;

    // a map of relative reproduction probabilities.
    string reproduction_file;

    /**
     * @brief Default constructor
     */
    SimParameters() : times(), configs()
    {
        fine_map_file = "none";
        coarse_map_file = "none";
        output_directory = "none";
        historical_fine_map_file = "none";
        historical_coarse_map_file = "none";
        sample_mask_file = "none";
        times_file = "null";
        dispersal_method = "none";
        landscape_type = "none";
        death_file = "none";
        reproduction_file = "none";
        dispersal_file = "none";
        min_speciation_gen = 0.0;
        max_speciation_gen = 0.0;
        is_protracted = false;
        restrict_self = false;
        m_prob = 0;
        cutoff = 0;
        tau = 0;
    }

    /**
     * @brief Links to the provided ConfigOption.
     * Assumes that the parameters have already been parsed from the config file.
     * @param configOption the parsed ConfigOption object
     */
    void importParameters(ConfigParser configOption)
    {
        configs = std::move(configOption);
        importParameters();
    }

    /**
     * @brief Imports the spatial variables from a path to the config file..
     * @param config_in string of the path to the config file
     */
    void importParameters(const string &conf_in)
    {
        // do the import of the values from combination of command-line arguments and file.
        configs.setConfig(conf_in, false, true);
        configs.parseConfig();
        importParameters();
    }

    /**
     * @brief Main import of parameters from the config file option.
     */
    void importParameters()
    {
        sample_x_size = stoul(configs.getSectionOptions("sample_grid", "x", "0"));
        sample_y_size = stoul(configs.getSectionOptions("sample_grid", "y", "0"));
        sample_x_offset = stoul(configs.getSectionOptions("sample_grid", "x_off", "0"));
        sample_y_offset = stoul(configs.getSectionOptions("sample_grid", "y_off", "0"));
        uses_spatial_sampling = static_cast<bool>(stoi(configs.getSectionOptions("sample_grid",
                                                                                 "uses_spatial_sampling", "0")));
        if(configs.hasSection("grid_map"))
        {
            grid_x_size = stoul(configs.getSectionOptions("grid_map", "x"));
            grid_y_size = stoul(configs.getSectionOptions("grid_map", "y"));
        }
        else
        {
            grid_x_size = sample_x_size;
            grid_y_size = sample_y_size;
        }
        sample_mask_file = configs.getSectionOptions("sample_grid", "path", "null");
        fine_map_file = configs.getSectionOptions("fine_map", "path", "none");
        fine_map_x_size = stoul(configs.getSectionOptions("fine_map", "x", "0"));
        fine_map_y_size = stoul(configs.getSectionOptions("fine_map", "y", "0"));
        fine_map_x_offset = stoul(configs.getSectionOptions("fine_map", "x_off", "0"));
        fine_map_y_offset = stoul(configs.getSectionOptions("fine_map", "y_off", "0"));
        coarse_map_file = configs.getSectionOptions("coarse_map", "path", "none");
        coarse_map_x_size = stoul(configs.getSectionOptions("coarse_map", "x", "0"));
        coarse_map_y_size = stoul(configs.getSectionOptions("coarse_map", "y", "0"));
        coarse_map_x_offset = stoul(configs.getSectionOptions("coarse_map", "x_off", "0"));
        coarse_map_y_offset = stoul(configs.getSectionOptions("coarse_map", "y_off", "0"));
        coarse_map_scale = stoul(configs.getSectionOptions("coarse_map", "scale", "0"));
        historical_fine_map_file = configs.getSectionOptions("historical_fine0", "path", "none");
        historical_coarse_map_file = configs.getSectionOptions("historical_coarse0", "path", "none");
        dispersal_method = configs.getSectionOptions("dispersal", "method", "none");
        m_prob = stod(configs.getSectionOptions("dispersal", "m_probability", "0"));
        cutoff = stod(configs.getSectionOptions("dispersal", "cutoff", "0.0"));
        // quick and dirty conversion for string to bool
        restrict_self = static_cast<bool>(stoi(configs.getSectionOptions("dispersal", "restrict_self", "0")));
        landscape_type = configs.getSectionOptions("dispersal", "landscape_type", "none");
        dispersal_file = configs.getSectionOptions("dispersal", "dispersal_file", "none");
        death_file = configs.getSectionOptions("death", "map", "none");
        reproduction_file = configs.getSectionOptions("reproduction", "map", "none");
        output_directory = configs.getSectionOptions("main", "output_directory", "Default");
        seed = stol(configs.getSectionOptions("main", "seed", "0"));
        job_type = stol(configs.getSectionOptions("main", "job_type", "0"));
        tau = stod(configs.getSectionOptions("main", "tau", "0.0"));
        sigma = stod(configs.getSectionOptions("main", "sigma", "0.0"));
        deme = stod(configs.getSectionOptions("main", "deme"));
        deme_sample = stod(configs.getSectionOptions("main", "sample_size"));
        max_time = stoul(configs.getSectionOptions("main", "max_time"));
        dispersal_relative_cost = stod(configs.getSectionOptions("main", "dispersal_relative_cost", "0"));
        spec = stod(configs.getSectionOptions("main", "min_spec_rate"));
        desired_specnum = stoul(configs.getSectionOptions("main", "min_species", "1"));
        if(configs.hasSection("protracted"))
        {
            is_protracted = static_cast<bool>(stoi(configs.getSectionOptions("protracted", "has_protracted", "0")));
            min_speciation_gen = stod(configs.getSectionOptions("protracted", "min_speciation_gen", "0.0"));
            max_speciation_gen = stod(configs.getSectionOptions("protracted", "max_speciation_gen"));
        }
        if(configs.hasSection("times"))
        {

            times_file = "set";
            auto times_str = configs.getSectionValues("times");
            for(auto i : times_str)
            {
                times.push_back(stod(i));
            }
            if(times.size() == 0)
            {
                times_file = "null";
            }
        }
        setHistorical(0);
    }

    /**
     * @brief Sets the main simulation parameters
     * @param job_type the job type reference number, used for file referencing
     * @param seed the seed to set random number generation
     * @param output_directory the output directory
     * @param max_time the maximum time to simulate for
     * @param desired_specnum the desired number of species to aim towards (currently not functional)
     * @param times_file the file containing a list of temporal sampling points
     */
    void setKeyParameters(const long long &job_type, const long long &seed, const string &output_directory,
                          const unsigned long &max_time, unsigned long desired_specnum,
                          const string &times_file)
    {
        this->job_type = job_type;
        this->seed = seed;
        this->output_directory = output_directory;
        this->max_time = max_time;
        this->desired_specnum = desired_specnum;
        this->times_file = times_file;

    }

    /**
     * @brief Sets the speciation parameters for the simulation.
     * @param spec_in the speciation rate to use
     * @param is_protracted_in if true, simulates as a protracted simulation
     * @param min_speciation_gen_in the minimum speciation generation for protracted simulations
     * @param max_speciation_gen_in the maximum speciation generation for protracted simulations
     */
    void setSpeciationParameters(const long double &spec_in, bool is_protracted_in, const double &min_speciation_gen_in,
                                 const double &max_speciation_gen_in)
    {
        spec = spec_in;
        is_protracted = is_protracted_in;
        min_speciation_gen = min_speciation_gen_in;
        max_speciation_gen = max_speciation_gen_in;
    }

    /**
     * @brief Sets the dispersal parameters for the simulation.
     * @param dispersal_method_in the method of individuals dispersing (normal, fat-tailed or norm-uniform)
     * @param sigma_in the sigma value for a normal distribution
     * @param tau_in the tau value for the fat-tailed distribution
     * @param m_prob_in the probability of uniform dispersal for the norm-uniform distribution
     * @param cutoff_in the maximum dispersal distance for the uniform distribution
     * @param dispersal_relative_cost_in the relative cost of dispersing through non-forest
     * @param restrict_self_in if true, prevents dispersal from the same cell
     * @param landscape_type_in the landscape type (infinite, tiled or closed)
     * @param dispersal_file_in a map of dispersal probabilities
     * @param reproduction_file_in a map of reproduction probabilities
     */
    void setDispersalParameters(const string &dispersal_method_in, const double &sigma_in, const double &tau_in,
                                const double &m_prob_in, const double &cutoff_in,
                                const double &dispersal_relative_cost_in, bool restrict_self_in,
                                const string &landscape_type_in, const string &dispersal_file_in,
                                const string &reproduction_file_in)
    {
        dispersal_method = dispersal_method_in;
        sigma = sigma_in;
        tau = tau_in;
        m_prob = m_prob_in;
        cutoff = cutoff_in;
        dispersal_relative_cost = dispersal_relative_cost_in;
        restrict_self = restrict_self_in;
        landscape_type = landscape_type_in;
        dispersal_file = dispersal_file_in;
        death_file = reproduction_file_in;
    }

    /**
     * @brief Sets the historical map parameters for the simulation.
     * @param historical_fine_file_map_in the fine resolution historical file
     * @param historical_coarse_map_file_in the coarse resolution historical file
     * @param gen_since_historical_in the number of generations since the historical state was achieved
     * @param habitat_change_rate_in the rate of habitat change towards the historical state
     */
    void setHistoricalMapParameters(const string &historical_fine_file_map_in,
                                    const string &historical_coarse_map_file_in,
                                    const double &gen_since_historical_in, const double &habitat_change_rate_in)
    {
        historical_fine_map_file = historical_fine_file_map_in;
        historical_coarse_map_file = historical_coarse_map_file_in;
        gen_since_historical = gen_since_historical_in;
        habitat_change_rate = habitat_change_rate_in;
    }

    void setHistoricalMapParameters(vector<string> path_fine, vector<unsigned long> number_fine,
                                    vector<double> rate_fine,
                                    vector<double> time_fine, vector<string> path_coarse,
                                    vector<unsigned long> number_coarse, vector<double> rate_coarse,
                                    vector<double> time_coarse)
    {
        habitat_change_rate = 0.0;
        if(!rate_fine.empty())
        {
            is_historical = true;
            habitat_change_rate = rate_fine[0];
        }
        gen_since_historical = 0.0;
        if(!time_fine.empty())
        {
            gen_since_historical = time_fine[0];
        }
        if(time_fine.size() != rate_fine.size() || rate_fine.size() != number_fine.size() ||
           number_fine.size() != path_fine.size())
        {
            stringstream ss;
            ss << "Lengths of historical fine map variables lists must be the same: " << time_fine.size() << "!=";
            ss << rate_fine.size() << "!=" << number_fine.size() << "!=" << path_fine.size() << endl;
            throw FatalException(ss.str());
        }
        if(time_coarse.size() != rate_coarse.size() || rate_coarse.size() != number_coarse.size() ||
           number_coarse.size() != path_coarse.size())
        {
            stringstream ss;
            ss << "Lengths of historical coarse map variables lists must be the same: " << time_coarse.size() << "!=";
            ss << rate_coarse.size() << "!=" << number_coarse.size() << "!=" << path_coarse.size() << endl;
            throw FatalException(ss.str());
        }
        for(unsigned long i = 0; i < time_fine.size(); i++)
        {
            string tmp = "historical_fine" + to_string(number_fine[i]);
            configs.setSectionOption(tmp, "path", path_fine[i]);
            configs.setSectionOption(tmp, "number", to_string(number_fine[i]));
            configs.setSectionOption(tmp, "time", to_string(time_fine[i]));
            configs.setSectionOption(tmp, "rate", to_string(rate_fine[i]));
        }
        for(unsigned long i = 0; i < time_coarse.size(); i++)
        {
            string tmp = "historical_coarse" + to_string(number_fine[i]);
            configs.setSectionOption(tmp, "path", path_coarse[i]);
            configs.setSectionOption(tmp, "number", to_string(number_coarse[i]));
            configs.setSectionOption(tmp, "time", to_string(time_coarse[i]));
            configs.setSectionOption(tmp, "rate", to_string(rate_coarse[i]));
        }
    }

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
                          const double &deme_in, const double &deme_sample_in, bool uses_spatial_sampling_in)
    {
        fine_map_file = fine_map_file_in;
        coarse_map_file = coarse_map_file_in;
        sample_mask_file = sample_mask_file_in;
        grid_x_size = grid_x_size_in;
        grid_y_size = grid_y_size_in;
        sample_x_size = sample_x_size_in;
        sample_y_size = sample_y_size_in;
        sample_x_offset = sample_x_offset_in;
        sample_y_offset = sample_y_offset_in;
        fine_map_x_size = fine_map_x_size_in;
        fine_map_y_size = fine_map_y_size_in;
        fine_map_x_offset = fine_map_x_offset_in;
        fine_map_y_offset = fine_map_y_offset_in;
        coarse_map_x_size = coarse_map_x_size_in;
        coarse_map_y_size = coarse_map_y_size_in;
        coarse_map_x_offset = coarse_map_x_offset_in;
        coarse_map_y_offset = coarse_map_y_offset_in;
        coarse_map_scale = coarse_map_scale_in;
        deme = deme_in;
        deme_sample = deme_sample_in;
        uses_spatial_sampling = uses_spatial_sampling_in;
    }

    /**
     * @brief Alters the historical parameters to the configuration matching the input number. If no configuration
     * option exists for this number, is_historical will be set to true.
     * @param n the historical map number to check.
     * @return bool true if we need to re-import the maps (i.e. the historical maps have changed between updates)
     */
    bool setHistorical(unsigned int n)
    {
        is_historical = true;
        bool finemapcheck = false;
        bool coarsemapcheck = false;
        // Loop over each element in the config file (each line) and check if it is historical fine or historical coarse.
        for(unsigned long i = 0; i < configs.getSectionOptionsSize(); i++)
        {
            if(configs[i].section.find("historical_fine") == 0)
            {
                // Then loop over each element to find the number, and check if it is equal to our input number.
                if(stol(configs[i].getOption("number")) == n)
                {
                    is_historical = false;
                    string tmpmapfile;
                    tmpmapfile = configs[i].getOption("path");
                    if(historical_fine_map_file != tmpmapfile)
                    {
                        finemapcheck = true;
                        historical_fine_map_file = tmpmapfile;
                    }
                    habitat_change_rate = stod(configs[i].getOption("rate"));
                    gen_since_historical = stod(configs[i].getOption("time"));
                }
            }
            else if(configs[i].section.find("historical_coarse") == 0)
            {
                if(stol(configs[i].getOption("number")) == n)
                {
                    string tmpmapfile;
                    tmpmapfile = configs[i].getOption("path");
                    is_historical = false;
                    if(tmpmapfile != historical_coarse_map_file)
                    {
                        coarsemapcheck = true;
                        historical_coarse_map_file = tmpmapfile;
                        // check matches
                        if(habitat_change_rate != stod(configs[i].getOption("rate")) ||
                           gen_since_historical != stod(configs[i].getOption("time")))
                        {
                            writeWarning(
                                    "Forest transform values do not match between fine and coarse maps. Using fine values.");
                        }
                    }
                }
            }
        }
        // if one of the maps has changed, we need to update, so return true.
        if(finemapcheck != coarsemapcheck)
        {
            return true;
        }
        else
        {
            // finemapcheck should therefore be the same as coarsemapcheck
            return finemapcheck;
        }
    }

    /**
     * @brief Prints selected important variables to the terminal.
     */
    void printVars()
    {
        stringstream os;
        os << "Seed: " << seed << endl;
        os << "Speciation rate: " << spec << endl;
        if(is_protracted)
        {
            os << "Protracted variables: " << min_speciation_gen << ", " << max_speciation_gen << endl;
        }
        os << "Job Type: " << job_type << endl;
        os << "Max time: " << max_time << endl;
        printSpatialVars();
        os << "Deme: " << deme << endl;
        os << "Deme sample: " << deme_sample << endl;
        os << "Output directory: " << output_directory << endl;
        os << "Disp Rel Cost: " << dispersal_relative_cost << endl;
        os << "Times: ";
        if(times_file == "set")
        {
            for(unsigned long i = 0; i < times.size(); i++)
            {
                os << times[i];
                if(i != times.size() - 1)
                {
                    os << ", ";
                }
            }
        }
        else
        {
            os << " 0.0";
        }
        os << endl;
        writeInfo(os.str());
    }

    /**
     * @brief Prints the spatial variables.
     */
    void printSpatialVars()
    {
        stringstream os;
        os << "Dispersal (tau, sigma): " << tau << ", " << sigma << endl;
        os << "Dispersal method: " << dispersal_method << endl;
        if(dispersal_method == "norm-uniform")
        {
            os << "Dispersal (m, cutoff): " << m_prob << ", " << cutoff << endl;
        }
        os << "Fine map\n-file: " << fine_map_file << endl;
        os << "-dimensions: (" << fine_map_x_size << ", " << fine_map_y_size << ")" << endl;
        os << "-offset: (" << fine_map_x_offset << ", " << fine_map_y_offset << ")" << endl;
        os << "Coarse map\n-file: " << coarse_map_file << endl;
        os << "-dimensions: (" << coarse_map_x_size << ", " << coarse_map_y_size << ")" << endl;
        os << "-offset: (" << coarse_map_x_offset << ", " << coarse_map_y_offset << ")" << endl;
        os << "-scale: " << coarse_map_scale << endl;
        os << "Sample grid" << endl;
        if(sample_mask_file != "none" && sample_mask_file != "null")
        {
            os << "-file: " << sample_mask_file << endl;
        }
        os << "-dimensions: (" << sample_x_size << ", " << sample_y_size << ")" << endl;
        os << "-optimised area: (" << grid_x_size << ", " << grid_y_size << ")" << endl;
        os << "-optimised offsets: (" << sample_x_offset << ", " << sample_y_offset << ")" << endl;
        writeInfo(os.str());
    }

    /**
     * @brief Sets the metacommunity parameters.
     * @param metacommunity_size the number of individuals in the community
     * @param speciation_rate the speciation rate for the metacommunity
     * @param seed the seed for the simulation
     * @param job_type the job referencing number
     */
    void setMetacommunityParameters(const unsigned long &metacommunity_size,
                                    const long double &speciation_rate,
                                    const unsigned long &seed,
                                    const unsigned long &job_type)
    {
        output_directory = "Default";
        // randomise the seed slightly so that we get a different starting number to the initial simulation
        this->seed = static_cast<long long int>(elegantPairing(seed, job_type));
        this->job_type = (long long int) job_type;
        deme = metacommunity_size;
        deme_sample = 1.0;
        spec = speciation_rate;
        // Default to 1000 seconds - should be enough for most simulation sizes, but can be changed later if needed.
        max_time = 1000;
        times_file = "null";
        min_speciation_gen = 0.0;
        max_speciation_gen = 0.0;
    }

    /**
     * @brief Overloading the << operator for outputting to the output stream
     * @param os the output stream.
     * @param m the SimParameters object.
     * @return os the output stream.
     */
    friend ostream &operator<<(ostream &os, const SimParameters &m)
    {
        os << m.fine_map_file << "\n" << m.coarse_map_file << "\n" << m.historical_fine_map_file << "\n";
        os << m.historical_coarse_map_file << "\n" << m.sample_mask_file << "\n";
        os << m.seed << "\n" << m.job_type << "\n" << m.grid_x_size << "\n" << m.grid_y_size << "\n";
        os << m.sample_x_size << "\n" << m.sample_y_size << "\n" << m.sample_x_offset << "\n" << m.sample_y_offset
           << "\n";
        os << m.fine_map_x_size << "\n" << m.fine_map_y_size << "\n";
        os << m.fine_map_x_offset << "\n" << m.fine_map_y_offset << "\n" << m.coarse_map_x_size << "\n"
           << m.coarse_map_y_size << "\n" << m.coarse_map_x_offset << "\n";
        os << m.coarse_map_y_offset << "\n" << m.coarse_map_scale << "\n" << m.desired_specnum << "\n";
        os << m.dispersal_relative_cost << "\n" << m.deme << "\n" << m.deme_sample << "\n";
        os << m.spec << "\n" << m.sigma << "\n" << m.max_time << "\n" << m.gen_since_historical << "\n"
           << m.habitat_change_rate << "\n" << m.tau;
        os << "\n" << m.dispersal_method << "\n";
        os << m.m_prob << "\n" << m.cutoff << "\n" << m.restrict_self << "\n" << m.landscape_type << "\n"
           << m.times_file << "\n";
        os << m.dispersal_file << "\n" << m.uses_spatial_sampling << "\n";
        os << m.times.size() << "\n";
        for(const auto &each : m.times)
        {
            os << each << "\n";
        }
        os << m.configs;
        return os;
    }

    /**
     * @brief Overloading the >> operator for inputting from an input stream
     * @param is the input stream
     * @param m the mapvars object
     * @return is the input stream
     */
    friend istream &operator>>(istream &is, SimParameters &m)
    {
        getline(is, m.fine_map_file);
        getline(is, m.coarse_map_file);
        getline(is, m.historical_fine_map_file);
        getline(is, m.historical_coarse_map_file);
        getline(is, m.sample_mask_file);
        is >> m.seed >> m.job_type >> m.grid_x_size >> m.grid_y_size;
        is >> m.sample_x_size >> m.sample_y_size >> m.sample_x_offset >> m.sample_y_offset;
        is >> m.fine_map_x_size >> m.fine_map_y_size;
        is >> m.fine_map_x_offset >> m.fine_map_y_offset >> m.coarse_map_x_size >> m.coarse_map_y_size
           >> m.coarse_map_x_offset;
        is >> m.coarse_map_y_offset >> m.coarse_map_scale >> m.desired_specnum >> m.dispersal_relative_cost >> m.deme
           >> m.deme_sample;
        is >> m.spec >> m.sigma >> m.max_time >> m.gen_since_historical >> m.habitat_change_rate >> m.tau;
        is.ignore();
        getline(is, m.dispersal_method);
        is >> m.m_prob >> m.cutoff >> m.restrict_self >> m.landscape_type;
        is.ignore();
        getline(is, m.times_file);
        getline(is, m.dispersal_file);
        is >> m.uses_spatial_sampling;
        unsigned long tmp_size;
        double tmp_time;
        is >> tmp_size;
        for(unsigned long i = 0; i < tmp_size; i++)
        {
            is >> tmp_time;
            m.times.push_back(tmp_time);
        }
        is >> m.configs;
        return is;
    }

    SimParameters &operator=(const SimParameters &other) = default;
};

#endif //SPECIATIONCOUNTER_SIMPARAMETERS_H
