# Example application of a metacommunity to a simulation
library(rcoalescence)


# First we run a simulation to generate a sample database

simulation <- TreeSimulation$new()
# Simulate 100 individuals
simulation$setSimulationParameters(seed=10, task =1, output_directory = "output", 
                                   min_speciation_rate = 0.5, deme=100)
simulation$runSimulation()
simulation$output()
