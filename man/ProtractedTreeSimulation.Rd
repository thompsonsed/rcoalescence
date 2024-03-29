% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppCoalescenceModule.R
\docType{class}
\name{ProtractedTreeSimulation}
\alias{ProtractedTreeSimulation}
\title{Protracted non-spatial neutral models}
\description{
Run non-spatial neutral models with protracted speciation. Here no map files are
required, but the protracted speciation parameters should be specified.
}
\section{Methods}{

\describe{
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

\section{Protracted non-spatial models}{

Run non-spatial neutral models with protracted speciation. Here no map files are required, but
the protracted speciation parameters should be specified.
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
}
