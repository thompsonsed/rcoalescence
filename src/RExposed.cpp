// This file is part of rcoalescence project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.


/**
 * @author Sam Thompson
 * @file RExposed.cpp
 * @brief Controls exposure of C++ objects to R using Rcpp.
 *
 * Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */

#include "Rcpp.h"
#include "RTreeSimulation.h"
#include "RSpatialTreeSimulation.h"
#include "RProtractedTreeSimulation.h"
#include "RProtractedSpatialTreeSimulation.h"
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(BH)]]
//' @useDynLib rcoalescence, .registration=TRUE

bool logging_mode;

RCPP_EXPOSED_CLASS(RTreeSimulation);

RCPP_EXPOSED_CLASS(RSpatialTreeSimulation);

RCPP_EXPOSED_CLASS(RProtractedTreeSimulation);

RCPP_EXPOSED_CLASS(RProtractedSpatialTreeSimulation);

RCPP_MODULE(coalescenceModule)
{
    using namespace Rcpp;
    class_<RTreeSimulation>("RTreeSimulation", "Simulates non-spatial neutral models.")
            .constructor("initialises the tree")
            .method("._setKeyParameters", &RTreeSimulation::setKeyParameters,
                    "Sets the main simulation parameters for house-keeping purposes.")
            .method("._setMinSpeciationRate", &RTreeSimulation::setMinSpeciationRate,
                    "Sets the minimum speciation rate for the simulation.")
            .method("._setSimulationProtractedParameters", &RTreeSimulation::setSimulationProtractedParameters,
                    "Sets the protracted variables for the simulation.")
            .method("._addMetacommunityParameters", &RTreeSimulation::addMetacommunityParameters,
                    "Adds the metacommunity parameters to the simulation.")
            .method("._addProtractedParameters", &RTreeSimulation::addProtractedParameters,
                    "Adds the protracted parameters for adding post-simulation.")
            .method("._addSpeciationRate", &RTreeSimulation::addSpeciationRate,
                    "Adds a speciaton rate for application post-simulation.")

            .method("setup", &RTreeSimulation::setup, "Completes all set-up routines for the simulation, "
                                                      "including assigning memory and importing the required files.")
            .method("runSimulation", &RTreeSimulation::runSimulation,
                    "Performs the actual simulation, returning true if the simulation is complete.")
            .method("._applySpeciationRates", &RTreeSimulation::applySpeciation,
                    "Applies the provided speciation parameters to the completed simulation.")
            .method("._getSpeciesAbundances", &RTreeSimulation::getSpeciesAbundances,
                    "Gets the species abundances for the last calculated set of speciation parameters.")
            .method("._getSpeciesRichness", &RTreeSimulation::getSpeciesRichness,
                    "Gets the species richness for the supplied community reference. ")
            .method("._getLastSpeciesRichness", &RTreeSimulation::getLastSpeciesRichness,
                    "Gets the last calculated species richness. ")
            .method("._output", &RTreeSimulation::output, "Writes the simulation to the output database.")
            .method("setLoggingMode", &RTreeSimulation::setLoggingMode,
                    "Sets the logging mode for the simulation object.")
            .method("getOutputDirectory", &RTreeSimulation::getOutputDirectory, "Gets the output directory.")
            .method("getSeed", &RTreeSimulation::getSeed, "Gets the simulation seed.")
            .method("getJobType", &RTreeSimulation::getJobType, "Gets the job type.")
            .method("._getSQLDatabase", &RTreeSimulation::getSQLDatabase, "Gets the SQL database path.")
            .method("._setDeme", &RTreeSimulation::setDeme, "Gets the deme size of the simulation.")
            .method("._setMultipleOutput", &RTreeSimulation::setMultipleOutput,
                    "Sets the multiple output variable, if it hasn't already been set.")
            .method("._getMultipleOutput", &RTreeSimulation::getMultipleOutput,
                    "Returns the value of the multiple_output variable.")
            .field("output_database", &RTreeSimulation::output_database);
    class_<RSpatialTreeSimulation>("RSpatialTreeSimulation", "Simulates spatially-explicit neutral models.")
            .constructor("initialises the spatial tree")
            .derives<RTreeSimulation>("RTreeSimulation")
            .method("._setDispersalParameters", &RSpatialTreeSimulation::setDispersalParameters,
                    "Sets the dispersal parameters for the simulation.")
            .method("._setHistoricalMapParameters", &RSpatialTreeSimulation::setHistoricalMapParameters,
                    "Sets the historical map parameters for the simulation.")
            .method("._setMapParameters", &RSpatialTreeSimulation::setMapParameters,
                    "Sets the map parameters for the simulation, including dimensions and offsets.")
            .method("._addHistoricalMap", &RSpatialTreeSimulation::addHistoricalMap,
                    "Adds a historical map with the desired parameter "
                    "set");
    class_<RProtractedTreeSimulation>("RProtractedTreeSimulation",
                                      "Simulates non-spatial neutral models with protracted speciation.")
            .constructor("initialises the tree")
            .derives<RTreeSimulation>("RTreeSimulation");
    class_<RProtractedSpatialTreeSimulation>("RProtractedSpatialTreeSimulation",
                                             "Simulates spatially-explicit neutral models with protracted speciation.")
            .constructor("initialises the spatial tree")
            .derives<RSpatialTreeSimulation>("RSpatialTreeSimulation");
}