# Example simulation using the supplied dimensions and offsets, with logging
library(rcoalescence)
simulation <- SpatialTree$new()
simulation$setSimulationParameters(seed=9, task =1, output_directory = "output", 
                                   speciation_rate = 0.5, sigma=2 * (2 ** 0.5), deme=1, 
                                   deme_sample = 0.1, 
                                   fine_map_file = "inst/extdata/sample/example_fine.tif",
                                   fine_map_x_size = 13, fine_map_y_size = 13, 
                                   coarse_map_file="inst/extdata/sample/example_coarse.tif",
                                   coarse_map_x_size = 35, coarse_map_y_size = 41, 
                                   coarse_map_x_offset = 11, coarse_map_y_offset = 14,
                                   sample_mask_file = "inst/extdata/sample/example_mask.tif", 
                                   sample_x_size = 13, sample_y_size = 13, uses_logging = TRUE)
simulation$runSimulation()

# Or rely on automatic detection of map dimensions and offsets, without logging
simulation <- SpatialTree$new()
simulation$setSimulationParameters(seed=9, task =1, output_directory = "output",
                                   speciation_rate = 0.001, 
                                   sigma=2 * (2 ** 0.5), deme=1, deme_sample = 0.1, 
                                   fine_map_file = "inst/extdata/sample/example_fine.tif", 
                                   coarse_map_file="inst/extdata/sample/example_coarse.tif",
                                   sample_mask_file = "inst/extdata/sample/example_mask.tif",
                                   uses_logging = FALSE)
simulation$runSimulation()

# Apply additional speciation rates
simulation$applySpeciationRates(speciation_rates = c(0.001, 0.7, 0.8))
# Output to a database
simulation$output()
# Get the species richness
simulation$getSpeciesRichness(1)
# Get the species abundances (last 5)
tail(simulation$getSpeciesAbundances(1), 5)


