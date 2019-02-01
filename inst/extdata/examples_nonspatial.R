# Example simulation for non-spatial models
library(rcoalescence)
# Create a new simulation object
simulation <- TreeSimulation$new()
# Simulate 100 individuals
simulation$setSimulationParameters(seed=10, task =1, output_directory = "output", 
                                   min_speciation_rate = 0.5, deme=100, uses_logging = FALSE)
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