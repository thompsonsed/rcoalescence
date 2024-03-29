% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppCoalescenceModule.R
\docType{class}
\name{ProtractedSpatialTreeSimulation}
\alias{ProtractedSpatialTreeSimulation}
\title{Protracted spatially-explicit neutral models}
\description{
Runs spatially-explicit neutral models with protracted speciation. Requires both
the protracted speciation parameters and map files defining spatial density.
}
\section{Methods}{

\describe{
\item{\code{applySpeciationRates(
  speciation_rates,
  output_file = "none",
  record_ages = FALSE,
  times_list = c(0),
  min_speciation_gens = c(0),
  max_speciation_gens = c(0),
  metacommunity_option = NA,
  metacommunity_size = NA,
  metacommunity_speciation_rate = NA,
  metacommunity_external_reference = NA
)}}{Applies the provided speciation parameters}

\item{\code{setSimulationParameters(
  task,
  seed,
  min_speciation_rate,
  output_directory = "output",
  max_time = 3600,
  times_list = c(0),
  uses_logging = NA,
  deme = 1,
  deme_sample = 1,
  min_speciation_gen = NA,
  max_speciation_gen = NA,
  partial_setup = FALSE
)}}{Sets all protracted simulation parameters}
}}

\section{Alternative classes}{

For alternative classes providing similar functionality, see \link[=NeutralTreeSimulation]{here}.
}

\section{Protracted spatial models}{

Runs spatially-explicit neutral models with protracted speciation. Requires both the protracted
speciation parameters and map files defining spatial density.
}

\section{Simulation parameters}{

These parameters are required at simulation time. The vast majority have defaults,
which are defined in setSimulationParameters()
\itemize{
\item \emph{task}: the task reference number
\item \emph{seed}: the seed for setting random number generation
\item \emph{output_directory}: the path to the folder for storing simulation files
\item \emph{max_time}: the maximum number of seconds to simulate for before pausing
\item \emph{times_list}: list of temporal sampling points
\item \emph{uses_logging}: if true, all outputs are written to console
\item \emph{deme}: the number of individuals per cell
\item \emph{deme_sample}: the global sampling proportion
\item \emph{min_speciation_rate}: the minimum speciation rate for the simulation
\item \emph{partial_setup}: if true, doesn't complete setup right away as additional steps must be performed. The user must call \code{setup()} manually.`
}
}

\section{Spatial parameters}{

These include dispersal parameters and maps files
\itemize{
\item \emph{sigma}: mean dispersal: distance for a normally-distributed dispersal kernel
\item \emph{dispersal_method}: the method of dispersal from "normal", "fat-tailed", "uniform" and "norm-uniform"
\item \emph{tau}: the tau parameter for a fat-tailed dispersal kernel
\item \emph{m_prob}: the probability of dispersing uniformly for the norm-uniform dispersal kernel
\item \emph{cutoff}: the maximum distance dispersal occurs for in a uniform dispersal kernel
\item \emph{dispersal_relative_cost}: the relative cost of moving through non-forest
\item \emph{restrict_self}: if true, prevents individuals from dispersing from their own cell
\item \emph{landscape_type}: type of landscape from "closed", "infinite" and "tiled"
\item \emph{dispersal_map_file}: a map of dispersal probabilities between pairs of cells in the landscape
\item \emph{reproduction_map_file}: a map of relative reproduction rates
\item \emph{death_map_file}: a map of relative death rates
\item \emph{fine_map_file}: fine resolution density map
\item \emph{coarse_map_file}: coarse resolution density map
\item \emph{sample_mask_file}: spatial sampling mask
\item \emph{grid_x_size}: x dimension of the grid
\item \emph{grid_y_size}: y dimension of the grid
\item \emph{sample_x_size}: x dimension of the sample mask
\item \emph{sample_y_size}: y dimension of the sample mask
\item \emph{sample_x_offset}: the x offset of the sample mask from the grid
\item \emph{sample_y_offset}: the y offset of the sample mask from the grid
\item \emph{fine_map_x_size}: the x dimension of the fine map
\item \emph{fine_map_y_size}: the y dimension of the fine map
\item \emph{fine_map_x_offset}: the x offset of the fine map from the sample mask
\item \emph{fine_map_y_offset}: the y offset of the fine map from the sample mask
\item \emph{coarse_map_x_size}: the x dimension of the coarse map
\item \emph{coarse_map_y_size}: the y dimension of the coarse map
\item \emph{coarse_map_x_offset}: the x offset of the coarse map from the fine map
\item \emph{coarse_map_y_offset}: the y offset of the coarse map from the fine map
\item \emph{coarse_map_scale}: the relative scale of the coarse map
\item \emph{uses_spatial_sampling}: if true, the sample mask defines relative sampling proportions across the map
\item \emph{historical_fine_map}: the historical fine map file
\item \emph{historical_coarse_map}: the historical coarse map file
\item \emph{gen_since_historical}: the number of generations since the historical state
\item \emph{habitat_change_rate}: the rate of change to the historical map
}
}

\section{Protracted speciation parameters}{

These parameters should match pairwise, i.e. the nth item in min_speciation_gens should go with
the nth item in max_speciation_gens
\itemize{
\item \emph{min_speciation_gens}: list of minimum number of generations required before speciation
\item \emph{max_speciation_gens}: list of maximum number of generations required before speciation
}
}

\section{Post-simulation parameters}{

These are for rebuilding the coalescence tree under different conditions.
\itemize{
\item \emph{output_file}: the directory to output to, defaults to "none"
\item \emph{use_spatial}: if true, records full spatial locations of all individuals. Default=FALSE
\item \emph{record_ages}: if true, records the ages of species to the database. Default=FALSE
\item \emph{sample_file}: supply a mask for defining spatial sampling. Default="null"
\item \emph{use_fragments}: supply a file containing fragment coordinates, or TRUE to let program calculate fragments
\item \emph{speciation_rates}: list of speciation rates to apply
\item \emph{times_list}: list of times to calculate coalescence tree at
\item \emph{min_speciation_gens}: list of the minimum number of generations required before speciation
\item \emph{max_speciation_gens}: list of the maximum number of generations required before speciation
}
}

\examples{
# Example simulation using the supplied dimensions and offsets, with logging
library(rcoalescence)
# Define the folder for the maps
# Define the directory containing the maps
getMapDir <- function() {
  return(system.file("sample", package = "rcoalescence"))
}
simulation <- ProtractedSpatialTreeSimulation$new()
simulation$setSimulationParameters(
  seed = 9,
  task = 1,
  output_directory = "output",
  min_speciation_rate = 0.5,
  sigma = 2 * (2**0.5),
  deme = 1,
  deme_sample = 0.1,
  fine_map_file = file.path(getMapDir(), "example_fine.tif"),
  fine_map_x_size = 13,
  fine_map_y_size = 13,
  coarse_map_file = file.path(getMapDir(), "example_coarse.tif"),
  coarse_map_x_size = 35,
  coarse_map_y_size = 41,
  coarse_map_x_offset = 11,
  coarse_map_y_offset = 14,
  sample_mask_file = file.path(getMapDir(), "example_mask.tif"),
  sample_x_size = 13,
  sample_y_size = 13,
  uses_logging = FALSE,
  min_speciation_gen = 100,
  max_speciation_gen = 1000
)
simulation$runSimulation()
# Or rely on automatic detection of map dimensions and offsets, without logging
simulation <- ProtractedSpatialTreeSimulation$new()
simulation$setSimulationParameters(
  seed = 9,
  task = 1,
  output_directory = "output",
  min_speciation_rate = 0.001,
  sigma = 2 * (2**0.5),
  deme = 1,
  deme_sample = 0.1,
  fine_map_file = file.path(getMapDir(), "example_fine.tif"),
  coarse_map_file = file.path(getMapDir(), "example_coarse.tif"),
  sample_mask_file = file.path(getMapDir(), "example_mask.tif"),
  uses_logging = FALSE,
  min_speciation_gen = 100,
  max_speciation_gen = 1000
)
simulation$runSimulation()
# Apply additional speciation rates
simulation$applySpeciationRates(
  speciation_rates = c(0.001, 0.7, 0.8),
  min_speciation_gens = c(10, 20, 20),
  max_speciation_gens = c(500, 500, 1000)
)
# Output to a database
simulation$output()
# Get the species richness
simulation$getSpeciesRichness(1)
# Get the species abundances (last 5)
tail(simulation$getSpeciesAbundances(1), 5)
# Remove the output directory
unlink("output", recursive = TRUE)
}
