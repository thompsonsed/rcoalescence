# Example simulation for non-spatial models with protracted speciation
library(rcoalescence)
simulation <- ProtractedTree$new()
# Simulate 100 individuals
simulation$setSimulationParameters(seed=11, task =1, output_directory = "output", 
                                   speciation_rate = 0.5, deme=100, uses_logging = TRUE,
                                   min_speciation_gen=100, max_speciation_gen = 1000)
simulation$runSimulation()

# Apply additional speciation rates
simulation$applySpeciationRates(speciation_rates = c(0.001, 0.7, 0.8),
                                min_speciation_gens = c(10, 20, 20), 
                                max_speciation_gens = c(500, 500, 1000))
# Output to a database
simulation$output()
# Get the species richness
simulation$getSpeciesRichness(1)
# Get the species abundances (last 5)
tail(simulation$getSpeciesAbundances(1), 5)

