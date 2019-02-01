# Example application of a metacommunity to a simulation
library(rcoalescence)


# First we run a simulation to generate a sample database
simulation <- TreeSimulation$new()
# Simulate 100 individuals in a non-spatial manner
simulation$setSimulationParameters(seed=10, task =1, output_directory = "output", 
                                   min_speciation_rate = 0.01, deme=100)
simulation$runSimulation()
# Apply some metacommunity parameters, each generating a separate community
simulation$setSpeciationParameters(speciation_rates = c(0.5, 0.9), 
                                   metacommunity_option = c("simulated", "analytical"), 
                                   metacommunity_size = c(10000, 1000000), 
                                   metacommunity_speciation_rate = c(0.001, 0.001))
simulation$applySpeciationRates()
# Output to a database (stored at output/data_10_1.db)
simulation$output()
# Print out the metacommunity parameters
simulation$getMetacommunityReferences()
simulation$getSpeciesRichness(1)
