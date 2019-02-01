# README #

*rcoalescence* is a sister project to [pycoalescence](http://pycoalescence.readthedocs.io/ "pycoalescence documentation") for running spatially-explicit neutral models for ecological systems.

### Key Features ###

* Simulate spatially-explicit neutral models with complex map configurations.
* Use varying spatial and temporal sampling regimes to mimic real data.
* Uses coalescence methods for high-performance simulations.
* Convenient, minimalistic code structure.
* Obtain a variety of biodiversity metrics including species richness, species' locations and 
  species abundances.

### Setup and installation ###

* First make sure the requirements are installed (see below).
* Then the package can be installed directly from bitbucket using ``devtools::install_bitbucket(thompsonsed/rcoalescence)``
* Currently also requires installation of ``necsim``, which is not yet publicly available - please 
contact authors for more information.

### Prerequisites ###

Before attempting installation, the following prequisites should be installed.

* The devtools package (for installing from bitbucket) and Rcpp (for compiling the package).
* GDAL, which can be installed through different methods, depending on your system. For macOS, use 
  ``brew install gdal``, for Ubuntu, ``apt-get install gdal-dev`` (or equivalent for other Linux
  distributions), and for Windows there are prebuilt binaries available
  [here](http://www.gisinternals.com/release.php).
* A C++ compiler for your system.
* The Boost library (also available with the [BH R package](https://cran.r-project.org/package=BH)).
* The Sqlite3 library, which usually comes with your system.

### Usage ###

WIP

### Contacts ###

* Sam Thompson (Imperial College London / National University Singapore)
	- thompsonsed@gmail.com
	- samuel.thompson14@imperial.ac.uk
* Based on ideas and code provided by James Rosindell (contact information on request).