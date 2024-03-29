# README #

[![CircleCI](https://img.shields.io/circleci/project/github/thompsonsed/rcoalescence.svg?label=CircleCI&logo=circleci)](https://circleci.com/github/thompsonsed/rcoalescence) | [![License](https://img.shields.io/badge/license-MIT-green)](https://opensource.org/licenses/MIT)


*rcoalescence* is a sister project to 
[pycoalescence](http://pycoalescence.readthedocs.io/ "pycoalescence documentation") for running spatially-explicit neutral models for ecological systems. The package is available for macOS and Linux operating systems.

*Package migrated from [bitbucket](https://bitbucket.org/thompsonsed/rcoalescence.git) on 07/03/2020.*

### Key Features ###

* Simulate spatially-explicit neutral models with complex map configurations.
* Use varying spatial and temporal sampling regimes to mimic real data.
* Uses coalescence methods for high-performance simulations.
* Convenient, minimalistic code structure.
* Obtain a variety of biodiversity metrics including species richness, species' locations and 
  species abundances.

### Setup and installation ###

#### Installing on Windows ####

- Make sure that you have installed R > 3.6.0 with [Rtools version 3.5 or later](https://cran.r-project.org/bin/windows/Rtools/index.html). 
- Install the necessary R package requirements using ``install.packages(c("devtools", "Rcpp", "RSQLite", "dplyr"))``.
- Run ``devtools::install_github("thompsonsed/rcoalescence")``. This will also download gdal (if it can't be found already).


#### Installing on macOS ####

- Ensure you have R > 3.4.0 installed.
- Install gdal >2.4.0. If you have homebrew installed, this can be done using ``brew install gdal2``. Check out the [brew site](https://brew.sh/) for the command to install brew.
- Make sure xcode is installed (``xcode-select --install``).
- Install the necessary R package requirements using ``install.packages(c("devtools", "Rcpp", "RSQLite", "dplyr"))``.
- Run ``devtools::install_github("thompsonsed/rcoalescence")``.

#### Installing on Linux ####

- Ensure you have R > 3.4.0 installed, plus a compiler that supports C++14 or later (e.g. gcc or clang).
- Install gdal >2.4.0 from your distribution's repositories (e.g. ``sudo apt-get install gdal-dev`` for Ubuntu), or build from source.
- Install the necessary R package requirements using ``install.packages(c("devtools", "Rcpp", "RSQLite", "dplyr"))``.
- Run ``devtools::install_github("thompsonsed/rcoalescence")``.

We also suggest a number of other packages contained within the **tidyverse** collection. All required packages can be installed using ``install.packages(c("devtools", "Rcpp", "RSQLite", "tidyverse", "rmarkdown", "knitr"))``.

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
  sample_mask_file = file.path(map_dir, "example_mask.tif"), # the path to the sample mask
  partial_setup=TRUE  # Don't fully import the maps yet - more to be added
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
