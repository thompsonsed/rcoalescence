// This file is part of rcoalescence project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.


/**
 * @author Sam Thompson
 * @file RExposed.cpp
 * @brief Controls logging from R.
 *
 * Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */

#include "Rcpp.h"
#include <Rcpp/r/headers.h>
#include "RTree.h"
#include "RSpatialTree.h"
#include "RProtractedTree.h"
#include "RProtractedSpatialTree.h"
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(BH)]]
//' @useDynLib rcoalescence, .registration=TRUE

bool logging_mode;


RCPP_EXPOSED_CLASS(RTree);
RCPP_EXPOSED_CLASS(RSpatialTree);
RCPP_MODULE(coalescenceModule) {
	using namespace Rcpp;
	class_<RTree>("RTree", "Simulates non-spatial neutral models.")
			.constructor("initialises the tree")
			.method("._setKeyParameters", &RTree::setKeyParameters,
					"Sets the main simulation parameters for house-keeping purposes.")
			.method("._setSpeciationParameters", &RTree::setSpeciationParameters,
					"Sets the speciation parameters for the simulation.")
			.method("setup", &RTree::setup, "Completes all set-up routines for the simulation, "
											"including assigning memory and importing the required files.")
			.method("runSimulation", &RTree::runSimulation,
					"Performs the actual simulation, returning true if the simulation is complete.")
			.method("._applySpeciationRates", &RTree::applySpeciation,
					"Applies the provided speciation parameters to the completed simulation.")
			.method("._getSpeciesAbundances", &RTree::getSpeciesAbundances,
					"Gets the species abundances for the last calculated set of speciation parameters.")
			.method("._getSpeciesRichness", &RTree::getSpeciesRichness, "Gets the species richness for the last "
																		"calculated set of speciation parameters.")
			.method("._output", &RTree::output, "Writes the simulation to the output database.")
			.method("setLoggingMode", &RTree::setLoggingMode, "Sets the logging mode for the simulation object.")
			.method("._getSQLDatabase", &RTree::getSQLDatabase)
			.method("._setDeme", &RTree::setDeme)
			.field("output_database", &RTree::output_database)
			;
	class_<RSpatialTree>("RSpatialTree", "Simulates spatially-explicit neutral models.")
			.constructor("initialises the spatial tree")
			.derives<RTree>("RTree")
			.method("._setDispersalParameters", &RSpatialTree::setDispersalParameters,
					"Sets the dispersal parameters for the simulation.")
			.method("._setHistoricalMapParameters", &RSpatialTree::setHistoricalMapParameters,
					"Sets the historical map parameters for the simulation.")
			.method("._setMapParameters", &RSpatialTree::setMapParameters,
					"Sets the map parameters for the simulation, including dimensions and offsets.")
			.method("._addHistoricalMap", &RSpatialTree::addHistoricalMap, "Adds a historical map with the desired parameter "
					"set")
			;
	class_<RProtractedTree>("RProtractedTree", "Simulates non-spatial neutral models with protracted speciation.")
			.constructor("initialises the tree")
			.derives<RTree>("RTree")
			;
	class_<RProtractedSpatialTree>("RProtractedSpatialTree", "Simulates spatially-explicit neutral models with protracted speciation.")
			.constructor("initialises the spatial tree")
			.derives<RProtractedTree>("RProtractedTree")
			.derives<RSpatialTree>("RSpatialTree")
			;
}
//
//RCPP_EXPOSED_CLASS(RTree);
//RCPP_MODULE(coalescenceModule) {
//	using namespace Rcpp;
//
//}
//
//RCPP_EXPOSED_CLASS(RProtractedSpatialTree);
//RCPP_MODULE(coalescenceModule) {
//	using namespace Rcpp;
//	class_<RProtractedSpatialTree>("RProtractedSpatialTree", "Simulates spatially-explicit neutral models with"
//														  " protracted speciation.")
//			.constructor("initialises the tree")
//			.method("._setKeyParameters", &RProtractedSpatialTree::setKeyParameters,
//					"Sets the main simulation parameters for house-keeping purposes.")
//			.method("._setSpeciationParameters", &RProtractedSpatialTree::setSpeciationParameters,
//					"Sets the speciation parameters for the simulation.")
//			.method("._setDispersalParameters", &RProtractedSpatialTree::setDispersalParameters,
//					"Sets the dispersal parameters for the simulation.")
//			.method("._setHistoricalMapParameters", &RProtractedSpatialTree::setHistoricalMapParameters,
//					"Sets the historical map parameters for the simulation.")
//			.method("._setMapParameters", &RProtractedSpatialTree::setMapParameters,
//					"Sets the map parameters for the simulation, including dimensions and offsets.")
//			.method("setup", &RProtractedSpatialTree::setup, "Completes all set-up routines for the simulation, "
//												   "including assigning memory and importing the required files.")
//			.method("._addHistoricalMap", &RProtractedSpatialTree::addHistoricalMap,
//					"Adds a historical map with the desired parameter set")
//			.method("runSimulation", &RProtractedSpatialTree::runSimulation, "Performs the actual simulation, "
//															"returning true if the simulation is complete.")
//			.method("._applySpeciationRates", &RProtractedSpatialTree::applySpeciation,
//					"Applies the provided speciation parameters to the completed simulation.")
//			.method("getSpeciesAbundances", &RProtractedSpatialTree::getSpeciesAbundances,
//					"Gets the species abundances for the last calculated set of speciation parameters.")
//			.method("._getSpeciesRichness", &RProtractedSpatialTree::getSpeciesRichness,
//					"Gets the species richness for the last calculated set of speciation parameters.")
//			.method("._output", &RProtractedSpatialTree::output, "Writes the simulation to the output database.")
//			.method("setLoggingMode", &RProtractedSpatialTree::setLoggingMode,
//					"Sets the logging mode for the simulation object.")
//			.method("._getSQLDatabase", &RProtractedSpatialTree::getSQLDatabase)
//			.field("output_database", &RProtractedSpatialTree::output_database)
//			;
//}
//
//
//RCPP_EXPOSED_CLASS(RProtractedTree);
//RCPP_MODULE(coalescenceModule) {
//	using namespace Rcpp;
//	class_<RProtractedTree>("RTree", "Simulates non-spatial neutral models with protracted speciation.")
//			.constructor("initialises the tree")
//			.method("._setKeyParameters", &RProtractedTree::setKeyParameters,
//					"Sets the main simulation parameters for house-keeping purposes.")
//			.method("._setSpeciationParameters", &RProtractedTree::setSpeciationParameters,
//					"Sets the speciation parameters for the simulation.")
//			.method("setup", &RProtractedTree::setup, "Completes all set-up routines for the simulation, "
//											"including assigning memory and importing the required files.")
//			.method("runSimulation", &RProtractedTree::runSimulation, "Performs the actual simulation, "
//															"returning true if the simulation is complete.")
//			.method("._applySpeciationRates", &RProtractedTree::applySpeciation, "Applies the provided speciation parameters"
//																	   "to the completed simulation.")
//			.method("getSpeciesAbundances", &RProtractedTree::getSpeciesAbundances, "Gets the species abundances for the "
//																		  "last calculated set of speciation parameters.")
//			.method("._getSpeciesRichness", &RProtractedTree::getSpeciesRichness, "Gets the species richness for the last "
//																		"calculated set of speciation parameters.")
//			.method("._output", &RProtractedTree::output, "Writes the simulation to the output database.")
//			.method("setLoggingMode", &RProtractedTree::setLoggingMode, "Sets the logging mode for the simulation object.")
//			.method("._getSQLDatabase", &RProtractedTree::getSQLDatabase)
//			.field("output_database", &RProtractedTree::output_database)
//			;
//}