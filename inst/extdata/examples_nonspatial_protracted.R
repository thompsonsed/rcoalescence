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