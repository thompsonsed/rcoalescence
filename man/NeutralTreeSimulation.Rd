% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppCoalescenceModule.R
\name{NeutralTreeSimulation}
\alias{NeutralTreeSimulation}
\title{Simulate neutral models on a variety of landscape conformations}
\description{
Four classes are provided here, \link{TreeSimulation}, \link{SpatialTreeSimulation},
\link{ProtractedTreeSimulation} and \link{ProtractedSpatialTreeSimulation}, for
non-spatial/spatial and non-protracted/protracted simulations. The process for setting up,
running and analysing simulations are the same, as described below.
}
\section{Non-spatial models}{

Run non-spatial neutral models. Here no map files are required, simply set up the
simulation and run. Dispersal parameters (sigma, tau, etc) are not relevant for this class.
}

\section{Spatial models}{

Run spatial neutral models. You must provide map files which provide the spatial density for the
simulation. Additional maps are optional.
}

\section{Protracted non-spatial models}{

Run non-spatial neutral models with protracted speciation. Here no map files are required, but
the protracted speciation parameters should be specified.
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

\section{Metacommunity parameters}{

These parameters are used for rebuilding the coalescence tree using a metacommunity, calculated
using either a non-spatial simulation, an analytical approximation, or an external database
which contains a community to draw individuals from.
\itemize{
\item \emph{metacommunity_option}: one of "simulated", "analytical" or a path to an external database
\item \emph{metacommunity_size}: the number of individuals in the metacommunity
\item \emph{metacommunity_speciation_rate}: the effective speciation rate of the metacommunity
\item \emph{external_reference}: an external reference for the Metacommunity.
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
# Example simulation for non-spatial models
library(rcoalescence)
# Create a new simulation object
simulation <- TreeSimulation$new()
# Simulate 100 individuals
simulation$setSimulationParameters(
  seed = 10, # the random number seed
  task = 1, # the job task for referencing files
  output_directory = "output", # the output directory
  min_speciation_rate = 0.5, # the minimum speciation rate used during simulation
  deme = 100, # the number of individuals to simulate
  uses_logging = FALSE # set to true to print out all simulation logs
)
# Run the actual simulation
simulation$runSimulation()
# Apply additional speciation rates
simulation$applySpeciationRates(speciation_rates = c(0.6, 0.7, 0.8))
# Output to a database
simulation$output()
# Get the species richness
simulation$getSpeciesRichness(1)
# Get the species abundances (last 5)
tail(simulation$getSpeciesAbundances(1), 5)
# Remove the output directory
unlink("output", recursive = TRUE)
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
# Example simulation for non-spatial models with protracted speciation
library(rcoalescence)
simulation <- ProtractedTreeSimulation$new()
# Simulate 100 individuals with protracted speciation occuring between generation 100 and 1000
# with a probability of 0.5 per generation
simulation$setSimulationParameters(
  seed = 11,
  task = 1,
  output_directory = "output",
  min_speciation_rate = 0.5,
  deme = 100,
  uses_logging = FALSE,
  min_speciation_gen = 100, # the minimum number of generations a lineage can exist for
  max_speciation_gen = 1000 # the maximum number of generations a lineage can exist for
)
simulation$runSimulation()
# Apply additional speciation rates
simulation$applySpeciationRates(
  speciation_rates = c(0.5, 0.7, 0.8),
  min_speciation_gens = c(10, 20, 20), # additional protracted minimums to apply
  max_speciation_gens = c(500, 500, 1000) # additional protracted maximums to apply
)
# Output to a database
simulation$output()
# Get the species richness
simulation$getSpeciesRichness(1)
# Get the species abundances (last 5)
tail(simulation$getSpeciesAbundances(1), 5)
# Remove the output directory
unlink("output", recursive = TRUE)
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
# Example application of a metacommunity to a simulation
library(rcoalescence)
# First we run a simulation to generate a sample database
simulation <- TreeSimulation$new()
# Simulate 100 individuals in a non-spatial manner
simulation$setSimulationParameters(
  seed = 10,
  task = 1,
  output_directory = "output",
  min_speciation_rate = 0.01,
  deme = 100
)
simulation$runSimulation()
# Apply some metacommunity parameters, each generating a separate community
simulation$applySpeciationRates(
  speciation_rates = c(0.5, 0.9), # the speciation rates to apply separately
  metacommunity_option = c("simulated", "analytical"), # the method for generating the metacommunity
  metacommunity_size = c(10000, 1000000), # the number of individuals in the metacommunity
  metacommunity_speciation_rate = c(0.001, 0.001) # the metacommunity speciation rate
)
# Output to a database (stored at output/data_10_1.db)
simulation$output()
# Get the community references (which have a column relating to the metacommunity parameters)
simulation$getCommunityReferences()
# Get the metacommunity references
simulation$getMetacommunityReferences()
# Get the species richness
simulation$getSpeciesRichness(1)
# Remove the output directory
unlink("output", recursive = TRUE)
}
