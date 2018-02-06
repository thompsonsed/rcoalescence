#' rcoalescence: Efficient spatially-explicit neutral ecology simulator written in C++
#' @include RcppExports.R
#' @importFrom Rcpp evalCpp
#' @name rcoalescence
#' @author Sam Thompson
#' @description A convenience wrapper for necsim, a Neutral Ecology Coalescence SIMulator for spatially-explicit neutral
#' models based in c++.
#' 
#' rcoalescence allows for large-scale simulations to be performed on fully spatial files, with features to suit a 
#' wide variety of use cases.
#' 
#' Methods are contained within the \linkS4class{SpatialTree} class for setting up, running and analysing
#' simulations. The general process is
#' \itemize{
#'     \item Set up the simulation with the provided SimulationParameters list and import required files.
#'     \item Run the simulation (generally the most time-consuming step).
#'     \item Apply speciation rates and other requirements post simulation using SpeciationParameters.
#'     \item Acquire the desired biodiversity metrics, such as species richness or species abundances.
#'     \item Optionally write results to an SQLite database for further manual analysis.
#'
#' }
#' 
#' @docType package
#' 
NULL

setRcppClass(Class="SpatialTree", CppClass="RSpatialTree", module="coalescenceModule", saveAs = "SpatialTree",
            methods = list(
              setKeyParameters = function(task, seed, output_directory="output",
                                          max_time=3600, desired_specnum=1,
                                          times_list=c(0.0)) {
                "Sets the main simulation parameters."
                ._setKeyParameters(task, seed, output_directory, max_time, desired_specnum, times_list)
                },
              setSpeciationParameters = function(speciation_rate, is_protracted=FALSE,
                                                 min_speciation_gen = 0.0, max_speciation_gen=0.0) {
                "Sets the speciation parameters for the simulation."
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
              setPristineMapParameters = function(pristine_fine_map="none", 
                                                  pristine_coarse_map="none", gen_since_pristine=100000000,
                                                  habitat_change_rate=0.0){
                "Sets the pristine map parameters for the simulation"
                ._setPristineMapParameters(pristine_fine_map, pristine_coarse_map, gen_since_pristine,
                                           habitat_change_rate)
                },
              setMapParameters = function(fine_map_file="null", coarse_map_file="none",
                               sample_mask_file="null", grid_x_size=10,
                               grid_y_size=10, sample_x_size=10,
                               sample_y_size=10, sample_x_offset=0,
                               sample_y_offset=0, fine_map_x_size=10,
                               fine_map_y_size=10,
                               fine_map_x_offset=0,
                               fine_map_y_offset=0,
                               coarse_map_x_size=10,
                               coarse_map_y_size=10,
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
                    grid_x_size <- grid_x_size
                    grid_y_size <- grid_y_size
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
              setParameters = function(p){
                "Sets all simulation parameters"
                setKeyParameters(p$task, p$seed, p$output_directory, p$max_time, p$desired_specnum, p$times_list)
                setSpeciationParameters(p$speciation_rate, p$is_protracted, p$min_speciation_gen, p$max_speciation_gen)
                setDispersalParameters(p$sigma, p$dispersal_method, p$tau, p$m_prob, p$cutoff, 
                                       p$dispersal_relative_cost, p$restrict_self, p$landscape_type,
                                       p$dispersal_file, p$reproduction_file)
                setMapParameters(p$fine_map_file, p$coarse_map_file,
                                 p$sample_mask_file, p$grid_x_size,
                                 p$grid_y_size, p$sample_x_size,
                                 p$sample_y_size, p$sample_x_offset,
                                 p$sample_y_offset, p$fine_map_x_size,
                                 p$fine_map_y_size,
                                 p$fine_map_x_offset,
                                 p$fine_map_y_offset,
                                 p$coarse_map_x_size,
                                 p$coarse_map_y_size,
                                 p$coarse_map_x_offset,
                                 p$coarse_map_y_offset,
                                 p$coarse_map_scale, p$deme,
                                 p$deme_sample, p$uses_spatial_sampling)
                setPristineMapParameters(p$pristine_fine_map, p$pristine_coarse_map, p$gen_since_pristine,
                                         p$habitat_change_rate)
                setup()},
              applySpeciationRates = function(p){
                "Applies the provided speciation parameters to the simulation"
                if(is.logical(p$use_fragments))
                {
                  p$use_fragments <- substr(as.character(p$use_fragments), 1, 1)
                }
                ._applySpeciationRates(p$output_file, p$use_spatial, p$sample_file,
                                  p$use_fragments, p$speciation_rates, p$times_list,
                                  p$min_speciation_gen, p$max_speciation_gen,
                                  p$metacommunity_size, p$metacommunity_speciation_rate)
                }
              )
            )
#' Container for speciation parameters, to be edited and supplied to the SpatialTree class 
#' during post-simulation coalescence tree building.
#' @name SpeciationParameters
#' @details Speciation parameters details.
#' * output_file the directory to output to, defaults to "none"
#' * use_spatial if true, records full spatial locations of all individuals. Default=FALSE
#' * sample_file supply a mask for defining spatial sampling. Default="null"
#' * use_fragments supply a file containing fragment coordinates, or TRUE to let program calculate fragments
#' * speciation_rates list of speciation rates to apply
#' * times_list list of times to calculate coalescence tree for
#' * min_speciation_gen the minimum number of generations required before speciation
#' * max_speciation_gen the maximum number of generations required before speciation
#' * metacommunity_size the number of individuals in the metacommunity
#' * metacommunity_speciation_rate the effective speciation rate of the metacommunity
#' @md
#' @export SpeciationParameters
SpeciationParameters <- list(output_file="none", use_spatial=FALSE, sample_file="null",
                             use_fragments=FALSE, speciation_rates=NA, times_list=c(0.0),
                             min_speciation_gen=0.0, max_speciation_gen=0.0,
                             metacommunity_size=0, metacommunity_speciation_rate=0.0)
#' Container for the simulation parameters, to be edited and supplied to the SpatialTree class for simulation.
#' @name SimulationParameteres
#' @details Simulation parameter details
#' * task the task reference number
#' * seed the seed for setting random number generation
#' * output_directory the path to the folder for storing simulation files
#' * max_time the maximum number of seconds to simulate for before pausing
#' * desired_specnum the desired number of species to aim for
#' * times_list list of temporal sampling points
#' * sigma mean dispersal distance for a normally-distributed dispersal kernel
#' * dispersal_method the method of dispersal from "normal", "fat-tailed", "uniform" and "norm-uniform"
#' * tau the tau parameter for a fat-tailed dispersal kernel
#' * m_prob the probability of dispersing uniformly for the norm-uniform dispersal kernel
#' * cutoff the maximum distance dispersal occurs for in a uniform dispersal kernel
#' * dispersal_relative_cost the relative cost of moving through non-forest
#' * restrict_self if true, prevents individuals from dispersing from their own cell
#' * landscape_type type of landscape from "closed", "infinite" and "tiled"
#' * dispersal_file a map of dispersal probabilities 
#' * reproduction_file a map of reproduction probabilies
#' * fine_map_file fine resolution density map
#' * coarse_map_file coarse resolution density map
#' * sample_mask_file spatial sampling mask
#' * grid_x_size x dimension of the grid
#' * grid_y_size y dimension of the grid
#' * sample_x_size x dimension of the sample mask
#' * sample_y_size y dimension of the sample mask
#' * sample_x_offset the x offset of the sample mask from the grid
#' * sample_y_offset the y offset of the sample mask from the grid
#' * fine_map_x_size the x dimension of the fine map
#' * fine_map_y_size the y dimension of the fine map
#' * fine_map_x_offset the x offset of the fine map from the sample mask
#' * fine_map_y_offset the y offset of the fine map from the sample mask
#' * coarse_map_x_size the x dimension of the coarse map
#' * coarse_map_y_size the y dimension of the coarse map
#' * coarse_map_x_offset the x offset of the coarse map from the fine map
#' * coarse_map_y_offset the y offset of the coarse map from the fine map
#' * coarse_map_scale the relative scale of the coarse map
#' * deme the number of individuals per cell
#' * deme_sample the global sampling proportion
#' * uses_spatial_sampling if true, the sample mask defines relative sampling proportions across the map
#' * pristine_fine_map the pristine fine map file
#' * pristine_coarse_map the pristine coarse map file
#' * gen_since_pristine the number of generations since the pristine state
#' * habitat_change_rate the rate of change to the pristine map
#' * speciation_rate the minimum speciation rate for the simulation
#' * is_protracted if true, simulation will be simulated with protracted speciation
#' * min_speciation_gen minimum number of generations required before speciation
#' * max_speciation_gen maximum number of generations required before speciation
#' @md
#' @export SimulationParameters
SimulationParameters <- list(task = NA, seed = NA, output_directory = "output",
                             max_time=3600, desired_specnum=1,
                             times_list=c(0.0), sigma = NA, dispersal_method = "normal",
                             tau=1.0, m_prob=0.0, cutoff=0,
                             dispersal_relative_cost=1.0, restrict_self=FALSE,
                             landscape_type="closed", dispersal_file="none",
                             reproduction_file="null",
                             fine_map_file="null", coarse_map_file="none",
                             sample_mask_file="null", grid_x_size=NA,
                             grid_y_size=NA, sample_x_size=10,
                             sample_y_size=10, sample_x_offset=0,
                             sample_y_offset=0, fine_map_x_size=10,
                             fine_map_y_size=10,
                             fine_map_x_offset=0,
                             fine_map_y_offset=0,
                             coarse_map_x_size=10,
                             coarse_map_y_size=10,
                             coarse_map_x_offset=0,
                             coarse_map_y_offset=0,
                             coarse_map_scale=1, deme=1,
                             deme_sample=1.0, uses_spatial_sampling=FALSE,
                             pristine_fine_map="none", 
                             pristine_coarse_map="none", gen_since_pristine=100000000,
                             habitat_change_rate=0.0, speciation_rate = NA, is_protracted=FALSE,
                             min_speciation_gen = 0.0, max_speciation_gen=0.0)


NULL