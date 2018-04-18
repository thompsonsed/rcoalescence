#' rcoalescence: Efficient spatially-explicit neutral ecology simulator written in C++
# @useDynLib rcoalescence, .registration=TRUE

#' @importFrom Rcpp evalCpp
#' @importFrom methods new
#' @import methods
#' @import RSQLite
#' @name rcoalescence
#' @author Sam Thompson
#' @description A convenience wrapper for necsim, a Neutral Ecology Coalescence SIMulator for 
#' spatially-explicit neutral models based in c++.
#' 
#' rcoalescence allows for large-scale simulations to be performed on fully spatial files, with 
#' features to suit a wide variety of use cases.
#' 
#' Methods are contained within the \link{SpatialTree} class for setting up, running and 
#' analysing simulations. The general process is
#' \itemize{
#'     \item Create a new SpatialTree object
#'     \item Set the simulation parameters and add additional historical maps.
#'     \item Run the simulation (generally the most time-consuming step).
#'     \item Apply speciation rates and other requirements post simulation using 
#'     SpeciationParameters.
#'     \item Get the species richness from the simulation.
#'     \item Optionally write results to an SQLite database for further manual analysis.
#'     \item If written to an output database, acquire the desired biodiversity metrics, such as 
#'     species locations or species abundances.
#'
#' }
#' 
#' @docType package
#' 
NULL

#' Simulate and analyse spatially-explicit neutral models.
#' @name SpatialTree
#' @export SpatialTree
#' @details Simulation parameters
#' @description Main class for setting up, running and analysing simulations.
#' These parameters are required at simulation time. The vast majority have defaults, 
#' which are defined in setSimulationParameters()
#' * *task*: the task reference number
#' * *seed*: the seed for setting random number generation
#' * *output_directory*: the path to the folder for storing simulation files
#' * *max_time*: the maximum number of seconds to simulate for before pausing
#' * *desired_specnum*: the desired number of species to aim for
#' * *times_list*: list of temporal sampling points
#' * *uses_logging*: if true, all outputs are written to console
#' * *sigma*: mean dispersal: distance for a normally-distributed dispersal kernel
#' * *dispersal_method*: the method of dispersal from "normal", "fat-tailed", "uniform" and "norm-uniform"
#' * *tau*: the tau parameter for a fat-tailed dispersal kernel
#' * *m_prob*: the probability of dispersing uniformly for the norm-uniform dispersal kernel
#' * *cutoff*: the maximum distance dispersal occurs for in a uniform dispersal kernel
#' * *dispersal_relative_cost*: the relative cost of moving through non-forest
#' * *restrict_self*: if true, prevents individuals from dispersing from their own cell
#' * *landscape_type*: type of landscape from "closed", "infinite" and "tiled"
#' * *dispersal_file*: a map of dispersal probabilities
#' * *reproduction_file*: a map of reproduction probabilies
#' * *fine_map_file*: fine resolution density map
#' * *coarse_map_file*: coarse resolution density map
#' * *sample_mask_file*: spatial sampling mask
#' * *grid_x_size*: x dimension of the grid
#' * *grid_y_size*: y dimension of the grid
#' * *sample_x_size*: x dimension of the sample mask
#' * *sample_y_size*: y dimension of the sample mask
#' * *sample_x_offset*: the x offset of the sample mask from the grid
#' * *sample_y_offset*: the y offset of the sample mask from the grid
#' * *fine_map_x_size*: the x dimension of the fine map
#' * *fine_map_y_size*: the y dimension of the fine map
#' * *fine_map_x_offset*: the x offset of the fine map from the sample mask
#' * *fine_map_y_offset*: the y offset of the fine map from the sample mask
#' * *coarse_map_x_size*: the x dimension of the coarse map
#' * *coarse_map_y_size*: the y dimension of the coarse map
#' * *coarse_map_x_offset*: the x offset of the coarse map from the fine map
#' * *coarse_map_y_offset*: the y offset of the coarse map from the fine map
#' * *coarse_map_scale*: the relative scale of the coarse map
#' * *deme*: the number of individuals per cell
#' * *deme_sample*: the global sampling proportion
#' * *uses_spatial_sampling*: if true, the sample mask defines relative sampling proportions across the map
#' * *historical_fine_map*: the historical fine map file
#' * *historical_coarse_map*: the historical coarse map file
#' * *gen_since_historical*: the number of generations since the historical state
#' * *habitat_change_rate*: the rate of change to the historical map
#' * *speciation_rate*: the minimum speciation rate for the simulation
#' * *is_protracted*: if true, simulation will be simulated with protracted speciation
#' * *min_speciation_gen*: minimum number of generations required before speciation
#' * *max_speciation_gen*: maximum number of generations required before speciation
#' @details Post-simulation parameter details
#' These are for rebuilding the coalescence tree under different conditions.
#' * *output_file*: the directory to output to, defaults to "none"
#' * *use_spatial*: if true, records full spatial locations of all individuals. Default=FALSE
#' * *sample_file*: supply a mask for defining spatial sampling. Default="null"
#' * *use_fragments*: supply a file containing fragment coordinates, or TRUE to let program calculate fragments
#' * *speciation_rates*: list of speciation rates to apply
#' * *times_list*: list of times to calculate coalescence tree for
#' * *min_speciation_gen*: the minimum number of generations required before speciation
#' * *max_speciation_gen*: the maximum number of generations required before speciation
#' * *metacommunity_size*: the number of individuals in the metacommunity
#' * *metacommunity_speciation_rate*: the effective speciation rate of the metacommunity
#' @md
#' @example inst/extdata/examples_1.R
#' 
SpatialTree <- setRcppClass("SpatialTree", "RSpatialTree", module="coalescenceModule",
                            fields=list(output_database = "character"),
             methods = list(
               setKeyParameters = function(task, seed, output_directory="output",
                                           max_time=3600, desired_specnum=1,
                                           times_list=c(0.0), uses_logging=NA) {
                 "Sets the key parameters for the simulation"
                 if(!is.na(uses_logging))
                 {
                   setLoggingMode(uses_logging)
                 }
                 ._setKeyParameters(task, seed, output_directory, max_time, desired_specnum, times_list)
                 },
               
               
               setSpeciationParameters = function(speciation_rate, is_protracted=FALSE,
                                                  min_speciation_gen = 0.0, max_speciation_gen=0.0) {
                 "Sets the speciation parameters for the simulation"
                 ._setSpeciationParameters(speciation_rate, is_protracted, min_speciation_gen, max_speciation_gen)
                 },
               
               
               setDispersalParameters = function(sigma, dispersal_method="normal", tau=1.0, m_prob=0.0, cutoff=0,
                                                 dispersal_relative_cost=1.0, restrict_self=FALSE,
                                                 landscape_type="closed", dispersal_file="null",
                                                 reproduction_file="null"){
                 "Sets the dispersal parameters for the simulation"
                 ._setDispersalParameters(sigma, dispersal_method, tau, m_prob, cutoff,
                                          dispersal_relative_cost, restrict_self, landscape_type,
                                          dispersal_file, reproduction_file)
                 },
               
               
               setMapParameters = function(fine_map_file="null", coarse_map_file="none",
                                           sample_mask_file="null", grid_x_size=NA,
                                           grid_y_size=NA, sample_x_size=0,
                                           sample_y_size=0, sample_x_offset=0,
                                           sample_y_offset=0, fine_map_x_size=0,
                                           fine_map_y_size=0,
                                           fine_map_x_offset=0,
                                           fine_map_y_offset=0,
                                           coarse_map_x_size=0,
                                           coarse_map_y_size=0,
                                           coarse_map_x_offset=0,
                                           coarse_map_y_offset=0,
                                           coarse_map_scale=1, deme=1,
                                           deme_sample=1.0, uses_spatial_sampling=FALSE){
                 "Sets the map parameters for the simulation"
                 if(sample_mask_file == "null")
                 {
                   sample_x_size <- fine_map_x_size
                   sample_y_size <- fine_map_y_size
                   if(is.na(grid_x_size) || is.na(grid_y_size))
                   {
                     grid_x_size <- fine_map_x_size
                     grid_y_size <- fine_map_y_size
                   }
                 }
                 else
                 {
                   if(is.na(grid_x_size) || is.na(grid_y_size))
                   {
                     grid_x_size <- sample_x_size
                     grid_y_size <- sample_y_size
                   }
                 }
                 ._setMapParameters(fine_map_file, coarse_map_file,
                                    sample_mask_file, grid_x_size,
                                    grid_y_size, sample_x_size,
                                    sample_y_size, sample_x_offset,
                                    sample_y_offset, fine_map_x_size,
                                    fine_map_y_size,
                                    fine_map_x_offset,
                                    fine_map_y_offset,
                                    coarse_map_x_size,
                                    coarse_map_y_size,
                                    coarse_map_x_offset,
                                    coarse_map_y_offset,
                                    coarse_map_scale, deme,
                                    deme_sample, uses_spatial_sampling)
               },
               
               
               setHistoricalMapParameters = function(historical_fine_map="none",
                                                   historical_coarse_map="none", gen_since_historical=100000000,
                                                   habitat_change_rate=0.0){
                 "Sets the historical map parameters for the simulation testing"
                 ._setHistoricalMapParameters(historical_fine_map, historical_coarse_map, gen_since_historical,
                                            habitat_change_rate)
                 },
               
               addHistoricalMap = function(historical_fine_map, historical_coarse_map="none", gen_since_historical=1,
                                         habitat_change_rate=0.0){
                 "Adds a historical map to the list of historical maps to use."
                 ._addHistoricalMap(historical_fine_map, historical_coarse_map, gen_since_historical,
                                  habitat_change_rate)
               },
               
               applySpeciationRates = function(speciation_rates, output_file="none", use_spatial=FALSE,
                                               sample_file="null", use_fragments=FALSE, times_list=c(0.0),
                                               min_speciation_gen=0.0, max_speciation_gen=0.0, 
                                               metacommunity_size=0, metacommunity_speciation_rate=0.0){
                 "Applies the provided speciation parameters to the simulation"
                 if(is.logical(use_fragments)){
                   use_fragments <- substr(as.character(use_fragments), 1, 1)
                   }
                 ._applySpeciationRates(output_file, use_spatial, sample_file,
                                        use_fragments, speciation_rates, times_list,
                                        min_speciation_gen, max_speciation_gen,
                                        metacommunity_size, metacommunity_speciation_rate)
               },
               
               
               setSimulationParameters = function(task, seed, speciation_rate, sigma,  output_directory="output",
                                                  max_time=3600, desired_specnum=1,
                                                  times_list=c(0.0), uses_logging=NA, is_protracted=FALSE,
                                                  min_speciation_gen = 0.0, max_speciation_gen=0.0,
                                                  dispersal_method="normal", tau=1.0, m_prob=0.0, cutoff=0,
                                                  dispersal_relative_cost=1.0, restrict_self=FALSE,
                                                  landscape_type="closed", dispersal_file="none",
                                                  reproduction_file="none", fine_map_file="null",
                                                  coarse_map_file="none",
                                                  sample_mask_file="null", grid_x_size=NA,
                                                  grid_y_size=NA, sample_x_size=0,
                                                  sample_y_size=0, sample_x_offset=0,
                                                  sample_y_offset=0, fine_map_x_size=0,
                                                  fine_map_y_size=0,
                                                  fine_map_x_offset=0,
                                                  fine_map_y_offset=0,
                                                  coarse_map_x_size=0,
                                                  coarse_map_y_size=0,
                                                  coarse_map_x_offset=0,
                                                  coarse_map_y_offset=0,
                                                  coarse_map_scale=1, deme=1,
                                                  deme_sample=1.0, uses_spatial_sampling=FALSE,
                                                  historical_fine_map="none",
                                                  historical_coarse_map="none", gen_since_historical=100000000,
                                                  habitat_change_rate=0.0
               ){
                 "Sets all simulation parameters"
                 setKeyParameters(task, seed, output_directory, max_time, desired_specnum, times_list, uses_logging)
                 setSpeciationParameters(speciation_rate, is_protracted, min_speciation_gen, max_speciation_gen)
                 setDispersalParameters(sigma, dispersal_method, tau, m_prob, cutoff,
                                        dispersal_relative_cost, restrict_self, landscape_type,
                                        dispersal_file, reproduction_file)
                 setMapParameters(fine_map_file, coarse_map_file,
                                  sample_mask_file, grid_x_size,
                                  grid_y_size, sample_x_size,
                                  sample_y_size, sample_x_offset,
                                  sample_y_offset, fine_map_x_size,
                                  fine_map_y_size,
                                  fine_map_x_offset,
                                  fine_map_y_offset,
                                  coarse_map_x_size,
                                  coarse_map_y_size,
                                  coarse_map_x_offset,
                                  coarse_map_y_offset,
                                  coarse_map_scale, deme,
                                  deme_sample, uses_spatial_sampling)
                 setHistoricalMapParameters(historical_fine_map, historical_coarse_map, gen_since_historical,
                                          habitat_change_rate)
                 setup()
               },
               
               
               setOutputDatabase = function(output_db){
                 "Checks the output database exists and has been created properly, if it has been set."
                 if(is.na(output_db) | output_db != "not_set"){
                   output_database <<- output_db
                   checkOutputDatabase()
                 }
                 else
                 {
                   stop("Output database has not been set.")
                 }
               },
               
               
               checkOutputDatabase = function(){
                 "Checks that the output database exists"
                 if(!file.exists(output_database)){
                   stop("Output database does not exist at " + output_database)
                 }
               },
               
               output = function(){
                 "Outputs the biodiversity data to an sql database."
                 ._output()
                 setOutputDatabase(._getSQLDatabase())
               },
               
               getCommunityReferences = function(){
                 "Returns a data frame containing all the community references and their parameter sets."
                 checkOutputDatabase()
                 conn <- dbConnect(SQLite(), output_database)
                 species_locations <- dbGetQuery(conn, "SELECT reference, speciation_rate, time, fragments, metacommunity_reference FROM COMMUNITY_PARAMETERS")
                 dbDisconnect(conn)
                 return(species_locations)
               },
               
               getMetacommunityReferences = function(){
                 "Returns a data frame containing all the metacommunity references and their parameter sets."
                 checkOutputDatabase()
                 conn <- dbConnect(SQLite(), output_database)
                 species_locations <- dbGetQuery(conn, "SELECT reference, speciation_rate, metacommunity_size FROM METACOMMUNITY_PARAMETERS")
                 dbDisconnect(conn)
                 return(species_locations)
               },
               
               getSpeciesLocations = function(community_reference = 1){
                 "Gets a data frame of species locations where the community reference matches the input."
                 checkOutputDatabase()
                 conn <- dbConnect(SQLite(), output_database)
                 species_locations <- dbGetQuery(conn, paste("SELECT species_id, x, y FROM SPECIES_LOCATIONS WHERE community_reference==", 
                                                             community_reference, sep=""))
                 dbDisconnect(conn)
                 return(species_locations)
               },
               
               getSpeciesAbundances = function(community_reference = 1){
                 "Gets a data frame of species abundances where the community reference matches the input"
                 checkOutputDatabase()
                 conn <- dbConnect(SQLite(), output_database)
                 species_locations <- dbGetQuery(conn, paste("SELECT species_id, no_individuals FROM SPECIES_ABUNDANCES WHERE community_reference==",
                                                             community_reference, sep=""))
                 dbDisconnect(conn)
                 return(species_locations)
               },
               
               getSpeciesRichness = function(community_reference=NA){
                 "Gets the community reference from the output database, or from the internal object if no 
                 community reference is supplied (this will return the last calculated species richness)."
                 if(is.na(community_reference)){
                   return(._getSpeciesRichness())
                 }
                 checkOutputDatabase()
                 conn <- dbConnect(SQLite(), output_database)
                 species_locations <- dbGetQuery(conn, paste("SELECT COUNT(DISTINCT(species_id)) FROM SPECIES_ABUNDANCES WHERE no_individuals > 0 AND ",
                                                             "community_reference ==", community_reference, sep=""))
                 dbDisconnect(conn)
                 return(species_locations[[1]])
               }
             )
             
)
# SpatialTree$methods(add=function(a){return(a+1)})
# print(RSpatialTree)
# RSpatialTree$add("add")
# spatial_gen$methods()