% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppCoalescenceModule.R
\docType{class}
\name{SpatialTreeSimulation}
\alias{SpatialTreeSimulation}
\title{Spatially-explicit neutral models}
\description{
Run spatial neutral models. You must provide map files which provide the spatial
density for the simulation. Additional maps are optional.
}
\section{Methods}{

\describe{
\item{\code{addHistoricalMap(
  historical_fine_map,
  historical_coarse_map = "none",
  gen_since_historical = 1,
  habitat_change_rate = 0
)}}{Adds a historical map to the list of historical maps to use.}

\item{\code{applySpeciationRates(
  speciation_rates = NA,
  output_file = "none",
  record_ages = FALSE,
  times_list = c(0),
  metacommunity_option = NA,
  metacommunity_size = NA,
  metacommunity_speciation_rate = NA,
  metacommunity_external_reference = NA
)}}{Applies the provided speciation parameters to the simulation.}

\item{\code{setDispersalParameters(
  sigma = 1,
  dispersal_method = "normal",
  tau = 1,
  m_prob = 0,
  cutoff = 0,
  dispersal_relative_cost = 1,
  restrict_self = FALSE,
  landscape_type = "closed",
  dispersal_map_file = "none"
)}}{Sets the dispersal parameters for the simulation}

\item{\code{setHistoricalMapParameters(
  historical_fine_map = "none",
  historical_coarse_map = "none",
  gen_since_historical = 1e+08,
  habitat_change_rate = 0
)}}{Sets the historical map parameters for the simulation.}

\item{\code{setMapParameters(
  fine_map_file = "null",
  coarse_map_file = "none",
  sample_mask_file = "null",
  grid_x_size = NA,
  grid_y_size = NA,
  sample_x_size = 0,
  sample_y_size = 0,
  sample_x_offset = 0,
  sample_y_offset = 0,
  fine_map_x_size = 0,
  fine_map_y_size = 0,
  fine_map_x_offset = 0,
  fine_map_y_offset = 0,
  coarse_map_x_size = 0,
  coarse_map_y_size = 0,
  coarse_map_x_offset = 0,
  coarse_map_y_offset = 0,
  coarse_map_scale = 1,
  deme = 1,
  deme_sample = 1,
  uses_spatial_sampling = FALSE,
  reproduction_map_file = "none",
  death_map_file = "none"
)}}{Sets the map parameters for the simulation}

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
  partial_setup = FALSE
)}}{Sets simulation parameters}
}}

\section{Alternative classes}{

For alternative classes providing similar functionality, see \link[=NeutralTreeSimulation]{here}.
}

\section{Spatial models}{

Run spatial neutral models. You must provide map files which provide the spatial density for the
simulation. Additional maps are optional.
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
# Example simulation using the supplied dimensions and offsets
library(rcoalescence)
# Define the directory containing the maps
getMapDir <- function() {
  return(system.file("sample", package = "rcoalescence"))
}

# The recommended method is to rely on rcoalescence to detect dimensions and offsets of your maps.
# Create a new object to contain the SpatialTreeSimulation
simulation <- SpatialTreeSimulation$new()
# Set the main simulation parameters
simulation$setSimulationParameters(
  seed = 1, # random number seed
  task = 10, # the job number for file naming
  output_directory = "output", # the output directory
  min_speciation_rate = 0.001, # the minimum speciation rate
  sigma = 2, # the dispersal sigma value
  deme = 10, # the number of individuals per cell (multiplied by the map values)
  deme_sample = 0.01, # the proportion of individuals sampled
  fine_map_file = file.path(getMapDir(), "example_fine.tif"), # the path to the fine resolution map
  coarse_map_file = file.path(getMapDir(), "example_coarse.tif"), # the path to the coarse map
  sample_mask_file = file.path(getMapDir(), "example_mask.tif"), # the path to the sample mask
  partial_setup=TRUE # Don't fully import files yet, as we will be adding more
)
# Add a historical map
simulation$addHistoricalMap(
  historical_fine_map = file.path(getMapDir(), "example_historical_fine.tif"),
  historical_coarse_map = file.path(getMapDir(), "example_coarse.tif")
)
# Run the actual simulation
simulation$runSimulation()

# Post-simulation
# Apply additional speciation rates
simulation$applySpeciationRates(speciation_rates = c(0.001, 0.7, 0.8))
# Output to a database
simulation$output()
# Get the species richness
simulation$getSpeciesRichness(1)
# Get the species abundances (last 5)
tail(simulation$getSpeciesAbundances(1), 5)

# Alternatively, specify all the dimensions, scaling and offsets of the maps manually.
simulation <- SpatialTreeSimulation$new()
simulation$setSimulationParameters(
  seed = 1, # random number seed
  task = 10, # the job number for file naming
  output_directory = "output", # the output directory
  min_speciation_rate = 0.001, # the minimum speciation rate
  sigma = 2, # the dispersal sigma value
  deme = 10, # the number of individuals per cell (multiplied by the map values)
  deme_sample = 0.01, # the proportion of individuals sampled
  fine_map_file = file.path(getMapDir(), "example_fine.tif"), # the path to the fine resolution map
  fine_map_x_size = 13, # the x dimension of the fine map
  fine_map_y_size = 13, # the y dimension of the fine map
  coarse_map_file = file.path(getMapDir(), "example_coarse.tif"), # the path to the coarse map
  coarse_map_x_size = 35, # x dimension of the coarse map
  coarse_map_y_size = 41, # the y dimension of the coarse map
  coarse_map_x_offset = 11, # the x offset of the coarse map
  coarse_map_y_offset = 14, # the y offset of the coarse map
  coarse_map_scale = 10.0, # the relative scale of the coarse map compared to the fine map
  sample_mask_file = file.path(getMapDir(), "example_mask.tif"), # the path to the sample mask
  sample_x_size = 13, # the x dimension of the sample file
  sample_y_size = 13 # the y dimension of the sample file
)

# Remove the output directory
unlink("output", recursive = TRUE)
}
