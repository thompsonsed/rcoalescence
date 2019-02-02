// This file is part of rcoalescence project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.
/**
 * @author Sam Thompson
 * @file RSpatialTree.cpp
 * @brief Wraps the SpatialTree class for running spatially-explicit neutral models from R.
 *
 * Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */
#include "RSpatialTreeSimulation.h"

RSpatialTreeSimulation::RSpatialTreeSimulation() : RTreeSimulation(), SpatialTree()
{
    output_database = "not_set";
    setLoggingMode(false);
}

RSpatialTreeSimulation::~RSpatialTreeSimulation()
= default;

void RSpatialTreeSimulation::setDispersalParameters(const double &sigma_in, const string &dispersal_method_in,
                                                    const double &tau_in, const double &m_prob_in,
                                                    const double &cutoff_in, const double &dispersal_relative_cost_in,
                                                    bool restrict_self_in, const string &landscape_type_in,
                                                    const string &dispersal_file_in, const string &reproduction_file_in)
{
    sim_parameters->setDispersalParameters(dispersal_method_in, sigma_in, tau_in, m_prob_in, cutoff_in,
                                           dispersal_relative_cost_in, restrict_self_in, landscape_type_in,
                                           dispersal_file_in, reproduction_file_in);
}

void RSpatialTreeSimulation::setHistoricalMapParameters(const string &historical_fine_file_map_in,
                                                        const string &historical_coarse_map_file_in,
                                                        const double &gen_since_historical_in,
                                                        const double &habitat_change_rate_in)
{
    sim_parameters->setHistoricalMapParameters(historical_fine_file_map_in, historical_coarse_map_file_in,
                                               gen_since_historical_in, habitat_change_rate_in);
}

void RSpatialTreeSimulation::addHistoricalMap(const string &fine_map, const string &coarse_map, const double &time,
                                              const double &rate)
{
    if(sim_parameters->historical_fine_map_file == "none")
    {
        sim_parameters->historical_fine_map_file = fine_map;
        sim_parameters->historical_coarse_map_file = coarse_map;
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

void RSpatialTreeSimulation::setMapParameters(const string &fine_map_file_in, const string &coarse_map_file_in,
                                              const string &sample_mask_file_in, const unsigned long &grid_x_size_in,
                                              const unsigned long &grid_y_size_in,
                                              const unsigned long &sample_x_size_in,
                                              const unsigned long &sample_y_size_in,
                                              const unsigned long &sample_x_offset_in,
                                              const unsigned long &sample_y_offset_in,
                                              const unsigned long &fine_map_x_size_in,
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
    sim_parameters->setMapParameters(fine_map_file_in, coarse_map_file_in, sample_mask_file_in, grid_x_size_in,
                                     grid_y_size_in, sample_x_size_in, sample_y_size_in, sample_x_offset_in,
                                     sample_y_offset_in, fine_map_x_size_in, fine_map_y_size_in, fine_map_x_offset_in,
                                     fine_map_y_offset_in, coarse_map_x_size_in, coarse_map_y_size_in,
                                     coarse_map_x_offset_in, coarse_map_y_offset_in, coarse_map_scale_in,
                                     deme_in, deme_sample_in, uses_spatial_sampling_in);
}

void RSpatialTreeSimulation::setup()
{
    if(!paths_fine.empty())
    {
        sim_parameters->setHistoricalMapParameters(paths_fine, numbers_fine, rates_fine, times_fine, paths_coarse,
                                                   numbers_coarse, rates_coarse, times_coarse);
    }
    SpatialTree::setup();
}

bool RSpatialTreeSimulation::runSimulation()
{
    if(!paths_fine.empty())
    {
        sim_parameters->setHistoricalMapParameters(paths_fine, numbers_fine, rates_fine, times_fine, paths_coarse,
                                                   numbers_coarse, rates_coarse, times_coarse);
        sim_parameters->setHistorical(0);
        landscape->resetHistorical();
    }
    return Tree::runSimulation();
}