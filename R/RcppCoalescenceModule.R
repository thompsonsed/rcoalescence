#' rcoalescence: Efficient spatial and non-spatial neutral ecology simulator written in C++
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
#' Four classes exist to represent non-spatial/spatial and non-protracted/protracted simulations. 
#' See \link[=NeutralTree]{here} for a full description.
#' 
#' @section Installation:
#' There are several dependencies required before installation of this package will be successful.
#' \itemize{
#'     \item Install boost library for c++, available \href{https://www.boost.org/}{here}.
#'     \item Install gdal library for c++, available \href{http://www.gdal.org/}{here}.
#'     \item Install the sqlite3 library for c++, available \href{https://www.sqlite.org/index.html}{here}.
#' }
#' 
#' @section Starting a simulation:
#' 
#' Methods are contained within the relevant classes \link[=NeutralTree]{(see here)} for setting up,
#' running and analysing simulations. The general process is
#' \itemize{
#'     \item Create a new \link{Tree}, \link{SpatialTree}, \link{ProtractedTree} or
#'           \link{ProtractedSpatialTree} object, depending on use case.
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
#' 
#' @docType package
#' 
NULL

#' Simulate neutral models on a variety of landscape conformations
#'
#' @name NeutralTree
#' @description Four classes are provided here, \link{Tree}, \link{SpatialTree}, 
#' \link{ProtractedTree} and \link{ProtractedSpatialTree}, for non-spatial/spatial and 
#' non-protracted/protracted simulations. The process for setting up, running and analysing 
#' simulations are the same, as described below.
#' 
#' @section Non-spatial Models:
#' Run non-spatial neutral models. Here no map files are required, simply set up the
#' simulation and run. Dispersal parameters (sigma, tau, etc) are not relevant for this class.
#' 
#' @section Spatial Models:
#' Run spatial neutral models. You must provide map files which provide the spatial density for the 
#' simulation. Additional maps are optional.
#' 
#' @section Protracted Non-spatial Models:
#' Run non-spatial neutral models with protracted speciation. Here no map files are required, but 
#' the protracted speciation parameters should be specified.
#' 
#' @section Protracted Spatial Models:
#' Runs spatially-explicit neutral models with protracted speciation. Requires both the protracted
#' speciation parameters and map files defining spatial density.
#' 
#' @section Simulation parameters:
#' These parameters are required at simulation time. The vast majority have defaults, 
#' which are defined in setSimulationParameters()
#' * *task*: the task reference number
#' * *seed*: the seed for setting random number generation
#' * *output_directory*: the path to the folder for storing simulation files
#' * *max_time*: the maximum number of seconds to simulate for before pausing
#' * *desired_specnum*: the desired number of species to aim for
#' * *times_list*: list of temporal sampling points
#' * *uses_logging*: if true, all outputs are written to console
#' * *deme*: the number of individuals per cell
#' * *deme_sample*: the global sampling proportion
#' * *speciation_rate*: the minimum speciation rate for the simulation
#' @section Spatial parameters:
#' These include dispersal parameters and maps files
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
#' * *uses_spatial_sampling*: if true, the sample mask defines relative sampling proportions across the map
#' * *historical_fine_map*: the historical fine map file
#' * *historical_coarse_map*: the historical coarse map file
#' * *gen_since_historical*: the number of generations since the historical state
#' * *habitat_change_rate*: the rate of change to the historical map
#' @section Protracted speciation parameters:
#' These parameters should match pairwise, i.e. the nth item in min_speciation_gens should go with
#' the nth item in max_speciation_gens
#' * *min_speciation_gens*: list of minimum number of generations required before speciation
#' * *max_speciation_gens*: list of maximum number of generations required before speciation
#' @section Post-simulation parameters:
#' These are for rebuilding the coalescence tree under different conditions.
#' * *output_file*: the directory to output to, defaults to "none"
#' * *use_spatial*: if true, records full spatial locations of all individuals. Default=FALSE
#' * *sample_file*: supply a mask for defining spatial sampling. Default="null"
#' * *use_fragments*: supply a file containing fragment coordinates, or TRUE to let program calculate fragments
#' * *speciation_rates*: list of speciation rates to apply
#' * *times_list*: list of times to calculate coalescence tree for
#' * *min_speciation_gens*: list of the minimum number of generations required before speciation
#' * *max_speciation_gens*: list of the maximum number of generations required before speciation
#' * *metacommunity_size*: the number of individuals in the metacommunity
#' * *metacommunity_speciation_rate*: the effective speciation rate of the metacommunity
#' @md
#' @example inst/extdata/examples_1.R
#' @example inst/extdata/examples_2.R
#' @example inst/extdata/examples_3.R
#' @example inst/extdata/examples_4.R
#' 
NULL

#' Non-spatial neutral simulations
#' @name Tree
#' @export Tree
#' @description Run non-spatial neutral models. Here no map files are required, simply set up the
#' simulation and run. Dispersal parameters (sigma, tau, etc) are not relevant for this class.
#' @section Alternative classes:
#' For alternative classes providing similar functionality, see \link[=NeutralTree]{here}.
#' @inheritSection NeutralTree Simulation parameters
#' @inheritSection NeutralTree Post-simulation parameters
Tree <- setRcppClass("Tree", "RTree", module="coalescenceModule",
                     fields=list(output_database = "character"),
                     methods = list(
                       setKeyParameters = function(task, seed, output_directory="output",
                                                   max_time=3600, desired_specnum=1,
                                                   times_list=c(0.0), uses_logging=NA, 
                                                   deme=1, deme_sample=1.0) {
                         "Sets the key parameters for the simulation"
                         if(!is.na(uses_logging))
                         {
                           setLoggingMode(uses_logging)
                         }
                         ._setKeyParameters(task, seed, output_directory, max_time, desired_specnum,
                                            times_list)
                         ._setDeme(deme, deme_sample)
                       },
                       
                       
                       setSpeciationParameters = function(speciation_rate) {
                         "Sets the speciation parameters for the simulation"
                         ._setSpeciationParameters(speciation_rate, FALSE, 0.0, 0.0)
                       },
                       
                       applySpeciationRates = function(speciation_rates, output_file="none", times_list=c(0.0),
                                                       metacommunity_size=0, metacommunity_speciation_rate=0.0){
                         "Applies the provided speciation parameters to the simulation"
                         ._applySpeciationRates(output_file, FALSE, "null",
                                                "FALSE", speciation_rates, times_list,
                                                c(0.0), c(0.0),
                                                metacommunity_size, metacommunity_speciation_rate)
                       },
                       
                       
                       setSimulationParameters = function(task, seed, speciation_rate, output_directory="output",
                                                          max_time=3600, desired_specnum=1,
                                                          times_list=c(0.0), uses_logging=NA,
                                                          deme=1,
                                                          deme_sample=1.0
                       ){
                         "Sets all simulation parameters"
                         setKeyParameters(task, seed, output_directory, max_time, 
                                          desired_specnum, times_list, uses_logging, deme, 
                                          deme_sample)
                         setSpeciationParameters(speciation_rate)
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
                         if(is.na(output_database)){
                           stop("Output database has not been set.")
                         }
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
                         species_locations = tryCatch(
                           {dbGetQuery(conn, 
                                       "SELECT reference, speciation_rate, time, fragments,
                                       metacommunity_reference, min_speciation_gen,
                                       max_speciation_gen FROM COMMUNITY_PARAMETERS")},
                           error= function(e) {
                             return(dbGetQuery(conn,
                                               "SELECT reference, speciation_rate, time, fragments,
                                               metacommunity_reference FROM COMMUNITY_PARAMETERS"))}
                         )
                         # species_locations <- dbGetQuery(conn, "SELECT reference, speciation_rate, time, fragments, metacommunity_reference FROM COMMUNITY_PARAMETERS")
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
                       
                       getSpeciesAbundances = function(community_reference = NA){
                         "Gets a data frame of species abundances where the community reference matches the input"
                         if(is.na(community_reference))
                         {
                           return(.getSpeciesAbundances())
                         }
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

#' Spatially-explicit neutral models
#' @name SpatialTree
#' @export SpatialTree
#' @description Run spatial neutral models. You must provide map files which provide the spatial 
#' density for the simulation. Additional maps are optional.
#' @inheritSection Tree Alternative classes
#' @inheritSection NeutralTree Spatial Models
#' @inheritSection NeutralTree Simulation parameters
#' @inheritSection NeutralTree Spatial parameters
#' @inheritSection NeutralTree Post-simulation parameters
#' @inheritParams Tree
SpatialTree <- setRcppClass("SpatialTree", "RSpatialTree", module="coalescenceModule",
                            fields=list(output_database = "character"),
                            contains=c("Tree"),
             methods = list(
               
               setDispersalParameters = function(sigma, dispersal_method="normal", tau=1.0,
                                                 m_prob=0.0, cutoff=0,
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
                                               metacommunity_size=0, metacommunity_speciation_rate=0.0){
                 "Applies the provided speciation parameters to the simulation"
                 if(is.logical(use_fragments)){
                   use_fragments <- substr(as.character(use_fragments), 1, 1)
                   }
                 ._applySpeciationRates(output_file, use_spatial, sample_file,
                                        use_fragments, speciation_rates, times_list,
                                        c(0.0), c(0.0),
                                        metacommunity_size, metacommunity_speciation_rate)
               },
               
               
               setSimulationParameters = function(task, seed, speciation_rate, sigma,  output_directory="output",
                                                  max_time=3600, desired_specnum=1,
                                                  times_list=c(0.0), uses_logging=NA, 
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
                 setKeyParameters(task, seed, output_directory, max_time, desired_specnum, 
                                  times_list, uses_logging)
                 setSpeciationParameters(speciation_rate)
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
               }
             )
             
)

#' Protracted non-spatial neutral models
#' @name ProtractedTree
#' @export ProtractedTree
#' @description Run non-spatial neutral models with protracted speciation. Here no map files are 
#' required, but the protracted speciation parameters should be specified.
#' @inheritSection Tree Alternative classes
#' @inheritSection NeutralTree Protracted Non-spatial Models
#' @inheritSection NeutralTree Simulation parameters
#' @inheritSection NeutralTree Protracted speciation parameters
#' @inheritSection NeutralTree Post-simulation parameters
#' @inheritParams Tree
ProtractedTree <- setRcppClass("ProtractedTree", "RProtractedTree", module="coalescenceModule",
                     fields=list(output_database = "character"),
                     contains=c("Tree"),
                     methods = list(
                       setSpeciationParameters = function(speciation_rate, 
                                                          min_speciation_gens = c(0.0),
                                                          max_speciation_gens=c(0.0)) {
                         "Sets the speciation parameters for the simulation"
                         ._setSpeciationParameters(speciation_rate, TRUE, min_speciation_gens,
                                                   max_speciation_gens)
                       },

                       applySpeciationRates = function(speciation_rates, output_file="none",
                                                       times_list=c(0.0),
                                                       min_speciation_gens=c(0.0), 
                                                       max_speciation_gens=c(0.0),
                                                       metacommunity_size=0, 
                                                       metacommunity_speciation_rate=0.0){
                         "Applies the provided speciation parameters to the simulation"
                         ._applySpeciationRates(output_file, FALSE, "null",
                                                "FALSE", speciation_rates, times_list,
                                                min_speciation_gens, max_speciation_gens,
                                                metacommunity_size, metacommunity_speciation_rate)
                       },


                       setSimulationParameters = function(task, seed, speciation_rate,
                                                          output_directory="output",
                                                          max_time=3600, desired_specnum=1,
                                                          times_list=c(0.0), uses_logging=NA,
                                                          min_speciation_gens = c(0.0),
                                                          max_speciation_gens = c(0.0),
                                                          deme=1,
                                                          deme_sample=1.0){
                         "Sets all simulation parameters"
                         setKeyParameters(task, seed, output_directory, max_time, desired_specnum,
                                          times_list, uses_logging, deme, deme_sample)
                         setSpeciationParameters(speciation_rate, min_speciation_gens,
                                                 max_speciation_gens)
                         setup()
                       }
                     )

)
#' Protracted spatially-explicit neutral models
#' @name ProtractedSpatialTree
#' @export ProtractedSpatialTree
#' @description  Runs spatially-explicit neutral models with protracted speciation. Requires both 
#' the protracted speciation parameters and map files defining spatial density.
#' @inheritSection Tree Alternative classes
#' @inheritSection NeutralTree Protracted Spatial Models
#' @inheritSection NeutralTree Simulation parameters
#' @inheritSection NeutralTree Spatial parameters
#' @inheritSection NeutralTree Protracted speciation parameters
#' @inheritSection NeutralTree Post-simulation parameters
#' @inheritParams Tree
#' @inheritParams SpatialTree
ProtractedSpatialTree <- setRcppClass("ProtractedSpatialTree", "RProtractedSpatialTree",
                                      module="coalescenceModule",
                               fields=list(output_database = "character"),
                               contains=c("ProtractedTree", "SpatialTree", "Tree"),
                               methods=list(
                                 setSimulationParameters = function(task, seed, speciation_rate, sigma,  output_directory="output",
                                                                    max_time=3600, desired_specnum=1,
                                                                    times_list=c(0.0), uses_logging=NA, 
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
                                                                    habitat_change_rate=0.0,
                                                                    min_speciation_gens=c(0.0),
                                                                    max_speciation_gens=c(0.0)
                                 ){
                                   "Sets all simulation parameters"
                                   setKeyParameters(task, seed, output_directory, max_time, desired_specnum, 
                                                    times_list, uses_logging)
                                   setSpeciationParameters(speciation_rate, min_speciation_gens=min_speciation_gens,
                                                           max_speciation_gens = max_speciation_gens)
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
                                 applySpeciationRates = function(speciation_rates, output_file="none", use_spatial=FALSE,
                                                                 sample_file="null", use_fragments=FALSE, times_list=c(0.0),
                                                                 metacommunity_size=0, metacommunity_speciation_rate=0.0, 
                                                                 min_speciation_gens=c(0.0), 
                                                                 max_speciation_gens=c(0.0)){
                                   "Applies the provided speciation parameters to the simulation"
                                   if(is.logical(use_fragments)){
                                     use_fragments <- substr(as.character(use_fragments), 1, 1)
                                   }
                                   ._applySpeciationRates(output_file, use_spatial, sample_file,
                                                          use_fragments, speciation_rates, times_list,
                                                          min_speciation_gens, max_speciation_gens,
                                                          metacommunity_size, metacommunity_speciation_rate)
                                 }
                               ))