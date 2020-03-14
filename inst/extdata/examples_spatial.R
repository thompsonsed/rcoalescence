# Example simulation using the supplied dimensions and offsets
library(rcoalescence)
# Define the directory containing the maps
getMapDir <- function(){
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
  sample_mask_file = file.path(getMapDir(), "example_mask.tif") # the path to the sample mask
)
# Add a historical map
simulation$addHistoricalMap(historical_fine_map =file.path(getMapDir(), "example_historical_fine.tif"),
                            historical_coarse_map = file.path(getMapDir(), "example_coarse.tif"))
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
