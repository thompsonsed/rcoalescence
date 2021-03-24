// This file is part of rcoalescence project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.
/**
 * @author Sam Thompson
 * @file RProtractedSpatialTreeSimulation.h
 * @brief Wraps the ProtractedSpatialTree class for running protracted spatially-explicit neutral models from R.
 *
 * Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */

#ifndef RCOALESCENCE_RPROTRACTEDSPATIALTREE_H
#define RCOALESCENCE_RPROTRACTEDSPATIALTREE_H

#ifndef WIN_INSTALL

#include <unistd.h>

#endif // WIN_INSTALL

#include <Rcpp.h>
#include "necsim/ProtractedSpatialTree.h"

using namespace necsim;

namespace rcoalescence
{
    /**
     * @brief Wrapper for spatial simulations.
     */
    class RProtractedSpatialTreeSimulation : public virtual RTreeSimulation, public virtual ProtractedSpatialTree
    {
    protected:
        vector<string> paths_fine;
        vector<unsigned long> numbers_fine;
        vector<double> rates_fine;
        vector<double> times_fine;
        vector<string> paths_coarse;
        vector<unsigned long> numbers_coarse;
        vector<double> rates_coarse;
        vector<double> times_coarse;
    public:
        RProtractedSpatialTreeSimulation() : RTreeSimulation(), ProtractedSpatialTree(), paths_fine(), numbers_fine(),
                                             rates_fine(), times_fine(), paths_coarse(), numbers_coarse(),
                                             rates_coarse(), times_coarse()
        {
            output_database = "not_set";
            setLoggingMode(false);
            bIsProtracted = true;

        }

        ~RProtractedSpatialTreeSimulation() override = default;

        using ProtractedSpatialTree::calcSpeciation;
        using ProtractedSpatialTree::speciateLineage;
        using ProtractedSpatialTree::getProtracted;
        using ProtractedSpatialTree::setProtractedVariables;
        using ProtractedSpatialTree::getProtractedVariables;
        using ProtractedSpatialTree::getProtractedGenerationMin;
        using ProtractedSpatialTree::getProtractedGenerationMax;
        using ProtractedSpatialTree::protractedVarsToString;

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
         */
        void setDispersalParameters(const double &sigma_in,
                                    const string &dispersal_method_in,
                                    const double &tau_in,
                                    const double &m_prob_in,
                                    const double &cutoff_in,
                                    const double &dispersal_relative_cost_in,
                                    bool restrict_self_in,
                                    const string &landscape_type_in,
                                    const string &dispersal_file_in)
        {
            sim_parameters->setDispersalParameters(dispersal_method_in,
                                                   sigma_in,
                                                   tau_in,
                                                   m_prob_in,
                                                   cutoff_in,
                                                   dispersal_relative_cost_in,
                                                   restrict_self_in,
                                                   landscape_type_in,
                                                   dispersal_file_in);
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
                                        const double &gen_since_historical_in,
                                        const double &habitat_change_rate_in)
        {
            sim_parameters->setHistoricalMapParameters(historical_fine_file_map_in,
                                                       historical_coarse_map_file_in,
                                                       gen_since_historical_in,
                                                       habitat_change_rate_in);
        }

        /**
         * @brief Adds the set of parameters for a historical map configuration.
         * @param fine_map the fine map file to use at the selected time
         * @param coarse_map the coarse map file to use at the selected time
         * @param time the time at which the map is relevant
         * @param rate the rate of change from the previous map to this map.
         */
        void addHistoricalMap(const string &fine_map, const string &coarse_map, const double &time, const double &rate)
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
         * @param reproduction_map_file_in the path to the map file of relative reproduction rates
         * @param death_map_file_in the path to the map file of relative death rates
         */
        void setMapParameters(const string &fine_map_file_in,
                              const string &coarse_map_file_in,
                              const string &sample_mask_file_in,
                              const unsigned long &grid_x_size_in,
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
                              const unsigned long &coarse_map_scale_in,
                              const unsigned long &deme_in,
                              const double &deme_sample_in,
                              bool uses_spatial_sampling_in,
                              const string &reproduction_map_file_in,
                              const string &death_map_file_in)
        {
            sim_parameters->setMapParameters(fine_map_file_in,
                                             coarse_map_file_in,
                                             sample_mask_file_in,
                                             grid_x_size_in,
                                             grid_y_size_in,
                                             sample_x_size_in,
                                             sample_y_size_in,
                                             sample_x_offset_in,
                                             sample_y_offset_in,
                                             fine_map_x_size_in,
                                             fine_map_y_size_in,
                                             fine_map_x_offset_in,
                                             fine_map_y_offset_in,
                                             coarse_map_x_size_in,
                                             coarse_map_y_size_in,
                                             coarse_map_x_offset_in,
                                             coarse_map_y_offset_in,
                                             coarse_map_scale_in,
                                             deme_in,
                                             deme_sample_in,
                                             uses_spatial_sampling_in);
            sim_parameters->setDemographicParameters(reproduction_map_file_in, death_map_file_in);
        }

        /**
         * @brief Calls SpatialTree::setup() to act as a wrapper accessible by R without extra classes.
         */
        void setupR() override
        {
            if(!paths_fine.empty())
            {
                sim_parameters->setHistoricalMapParameters(paths_fine,
                                                           numbers_fine,
                                                           rates_fine,
                                                           times_fine,
                                                           paths_coarse,
                                                           numbers_coarse,
                                                           rates_coarse,
                                                           times_coarse);
            }
            SpatialTree::setup();
        }

        /**
         * @brief Calls SpatialTree::runSimulation() to act as a wrapper accessible by R without extra classes.
         * @return bool true if simulation completes successfully
         */
        bool runSimulation()
        {
            if(!paths_fine.empty())
            {
                sim_parameters->setHistoricalMapParameters(paths_fine,
                                                           numbers_fine,
                                                           rates_fine,
                                                           times_fine,
                                                           paths_coarse,
                                                           numbers_coarse,
                                                           rates_coarse,
                                                           times_coarse);
                sim_parameters->setHistorical(0);
                landscape->resetHistorical();
            }
            return Tree::runSimulation();
        }
    };
}
#endif //RCOALESCENCE_RPROTRACTEDSPATIALTREE_H
