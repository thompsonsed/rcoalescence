# README #

[![CircleCI](https://circleci.com/bb/thompsonsed/rcoalescence.svg?style=svg)](https://circleci.com/bb/thompsonsed/rcoalescence)

*rcoalescence* is a sister project to 
[pycoalescence](http://pycoalescence.readthedocs.io/ "pycoalescence documentation") for running spatially-explicit neutral models for ecological systems. The package is available for macOS and Linux operating systems.

### Key Features ###

* Simulate spatially-explicit neutral models with complex map configurations.
* Use varying spatial and temporal sampling regimes to mimic real data.
* Uses coalescence methods for high-performance simulations.
* Convenient, minimalistic code structure.
* Obtain a variety of biodiversity metrics including species richness, species' locations and 
  species abundances.

### Setup and installation ###

* First make sure the requirements are installed (see below).
* Then the package can be installed directly from bitbucket using ``devtools::install_bitbucket("thompsonsed/rcoalescence")``
* Currently also requires installation of ``necsim``, which is not yet publicly available - please 
contact authors for more information.

### Prerequisites ###

Before attempting installation, the following prequisites should be installed.

* The devtools package (for installing from bitbucket) and Rcpp (for compiling the package). Run ``install.packages(c("devtools", "Rcpp"))``
* GDAL, which can be installed through different methods, depending on your system. For macOS, use ``brew install gdal``, for Ubuntu, ``apt-get install gdal-dev`` (or equivalent for other Linux distributions), and for Windows there are prebuilt binaries available [here](http://www.gisinternals.com/release.php).
* A C++ compiler for your system. On macOS, make sure xcode is installed (``xcode-select --install``). Linux systems should come with gcc by default.
* The Boost library (also available with the [BH R package](https://cran.r-project.org/package=BH)).
* The Sqlite3 library, which usually comes with your system.

### Usage ###

See the examples by running ``?rcoalescence`` after importing the package into R. The basic process 
for a spatial simulation is:

```R
library(rcoalescence)

# Define the folder containing the maps (the examples come with the package)
map_dir <- system.file("sample", package = "rcoalescence")

# The recommended method is to rely on rcoalescence to detect dimensions and offsets of your maps.
# For specifying all dimensions and offsets manually, please see the examples.
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
  deme_sample = 0.1, # the proportion of individuals sampled
  fine_map_file = file.path(map_dir, "example_fine.tif"), # the path to the fine resolution map
  coarse_map_file = file.path(map_dir, "example_coarse.tif"), # the path to the coarse  map
  sample_mask_file = file.path(map_dir, "example_mask.tif") # the path to the sample mask
)
# Add a historical map
simulation$addHistoricalMap(historical_fine_map =file.path(map_dir, "example_historical_fine.tif"),
                            historical_coarse_map = file.path(map_dir, "example_coarse.tif"))
# Run the actual simulation
simulation$runSimulation() 

# Post-simulation
# Apply additional speciation rates
simulation$applySpeciationRates(speciation_rates = c(0.001, 0.7, 0.8))
# Output to a database
simulation$output()
# Get the species richness
simulation$getSpeciesRichness(1)
```

### Contacts ###

* Sam Thompson (Imperial College London / National University Singapore)
	- thompsonsed@gmail.com
	- samuel.thompson14@imperial.ac.uk
* Based on ideas and code provided by James Rosindell (contact information on request).
