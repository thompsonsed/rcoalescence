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
