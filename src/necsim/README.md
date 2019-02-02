# necsim #

Version: 1.0
This project is released under MIT licence
See file **LICENSE.txt** or go to [here](https://opensource.org/licenses/MIT) for full license details.

## CONTENTS ##
* **INTRODUCTION**
* **INSTRUCTIONS**
* **REQUIREMENTS**
* **DEGUGGING**
* **CLASS DESCRIPTIONS**
* **KNOWN BUGS**
* **FAQS**
* **CONTACTS**

## INTRODUCTION ##

necsim is a generic spatial coalescence simulator for neutral systems. It applies the model to temporally and spatially varying density maps for a specific set of supplied current_metacommunity_parameters, and outputs information for each individual to a SQL database.

necsim includes functionality for applying varying speciation rates after simulations are complete. This enables the main simulation to be run with the *minimum* speciation rate required and afterwards analysis can be completed using different speciation rates.

The recommended method of usage is through the pycoalescence package, using a Python interface for installation, simulation setup and running. See [here](http://pycoalescence.readthedocs.io/) for more details.

You are free to modify and distribute the code as per the license specified in **LICENCE.txt** to suit any additional neutral simulation requirements (or any other purpose).

## INSTRUCTIONS ##
###Compiling the program###
See the Requirements section for a full species_id_list of the necessary prerequisites. Once these are installed, compiling the program should be relatively easy. necsim requires a linker to the boost libraries, as well as the sqlite3 library. It is recommended to run with the maximum optimisation possible.


Additionally, if support is required for tif files (an alternative to importing csv files), the [gdal library](http://www.gdal.org/) is required. See the online documentation for help compiling gdal for your operating system. When compiling using gdal, use the ```-D with_gdal``` compilation flag.

For compilation on High Performance Computing (HPC) systems, they will likely use intel compilers. The header files for the sqlite and boost packages may need to be copied in to the working directory to avoid problems with linking to libraries. Check the service providers' documentation for whether these libraries are already installed on the HPC. 
for the application of different speciation rates.

###Running simulations###
Note that the recommended method of running and analysing simulations is through the [**pycoalescence**](https://pycoalescence.readthedocs.io) Python package.
The routine relies on supplying command line arguments (see below) for all the major simulation variables. Alternatively, supplying a config .txt file and using the command line arguments `./necsim -c /path/to/config.txt` can be used for parsing command line arguments from the text file. 

####Command Line Arguments ####
The following command line arguments are required. This species_id_list can be accessed by running `“./necsim -h”` or `./necsim -help`

The command line options to be specified are:

1. the seed for the simulation.
2. the simulation job_type (for file reference).
3. the map config file.
4. the output directory.
5. the minimum speciation rate.
6. the dispersal sigma value.
7. the dispersal tau value.
8. the deme size.
9. the deme sample size.
10. the maximum simulation time (in seconds).
11. the lambda value for moving through non-habitat.
12. the temporal sampling file containing tab-separated generation values for sampling points in time (null for only sampling the present)
13. the minimum number of species known to exist. (Currently has no effect).
14. (and onwards) speciation rates to apply after simulation.

In this set up, the map config file contains a file on each line, with tab separation between the different variables. The "ref" flag contains the object type, followed by all other current_metacommunity_parameters. An example is given below.

ref=sample_grid	path=/path/to/file	x=100	y=200	mask=/path/to/mask
ref=fine_map	path=/path/to/file	x=100	y=200	x_off=10	y_off=20
ref=historical_fine	path=/path/to/file	number=n	rate=r	time=g

Alternatively, by specifying the -f flag, (full mode) as the first argument, the program can read in extended command line arguments, which are as followed.

1. the number used for setting the seed.
2. the sample grid x dimension
3. the sample grid y dimension
4. the fine map file relative path.
5. the fine map x dimension
6. the fine map y dimension
7. the fine map x offset
8. the fine map y offset
9. the coarse map file relative path.
10. the coarse map x dimension
11. the coarse map y dimension
12. the coarse map x offset
13. the coarse map y offset
14. the scale of the coarse map compared to the fine (10 means resolution of coarse map = 10 x resolution of fine map)
15. the output directory
16. the speciation rate.
17. the dispersal sigma value.
18. the deme size
19. the deme sample size (as a proportion of deme size)
20. the time to run the simulation (in seconds).
21. lambda - the relative cost of moving through non-forest
22. job_type - for referencing the specific job_type later on.
23. the minimum number of species the system is known to contain.
24. the historical fine map file to use
25. the historical coarse map file to use
26. the rate of forest change from historical
27. the time (in generations) since the historical forest was seen.
28. the dispersal tau value (the width of the kernel.
29. the sample mask, with binary 1:0 values for areas that we want to sample from. If this is not provided then this will default to mapping the whole area.
30.  the link to the file containing every generation that the species_id_list should be expanded. This should be in the format of a species_id_list.
31. (and onwards) - speciation rates to apply after the simulation is complete.

####Config Files ####
The program also accepts a config file, specified by running `./necsim -c /path/to/config.txt`. The format of the config file is
```
rand_seed = i
sample_x_dim = i
sample_y_dim = i
fine_source = /path/to/fine.csv
fine_x_dim = i
fine_y_dim = i
fine_x_offset = i
fine_y_offset = i
coarse_source = /path/to/coarse.csv
coarse_x_dim = i
coarse_y_dim = i
coarse_x_offset = i
coarse_y_offset = i
coarse_scale = i
output_dir = /path/to/outdir
spec_rate = d
zfat = f
deme_size = i
deme_sample = d
wall_time = i
lambda = 1
job_num = i
est_spec = i
historical_fine_source = /path/to/historical/fine.csv
historical_coarse_source = /path/to/historical/coarse.csv
forest_change = d
time_since = f
dispersal = f
sampledatamask = /path/to/sample/mask.csv
time_config_file = /path/to/time/file.txt
speciationrate1 = d
speciationrate2 = d
...
```
where `i` represents a positive integer, `d` is a decimal value between 0 and 1, and `f` is any positive number (float). Whilst this does help with readability of the code, the order of the arguments is essential at this stage (i.e. don't switch the order of the lines). Future versions may alter the system of reading such that the current_metacommunity_parameters are set according to their key. Any number of speciation rates (or 0) can be at the end of the file.

####Outputs####
Upon successful completion of a simulation, necsim will produce an SQLite database file in the output directory in an SQL\_data folder. This database contains several tables, which can be accessed using a program like [DB Browser for SQLite](http://sqlitebrowser.org/) or Microsoft Access. Alternatively, most programming languages have an SQLite interface ([RSQlite](https://cran.r-project.org/web/packages/RSQLite/index.html), [Python sqlite3](https://docs.python.org/2/library/sqlite3.html))

* The main table within the database is the SPECIES\_LIST table, which is the location and inheritence of every lineage recorded. Several other important data structures (such as whether it is a "tip" of the phylogenetic tree of not) which are used  when re-constructing the species identity.
* A secondary output from necsims is a SIMULATION\_PARAMETERS table for identifying the exact current_metacommunity_parameters with which the model is run.
* SpeciationCounter also produces a SPECIES_ABUNDANCES table containing species abundances across the whole sample map, plus (optionally) a table of SPECIES\_LOCATIONS (containing the x,y location of every individual) and FRAGMENT\_ABUNDANCES (species abundances for each habitat fragment separately).

## REQUIREMENTS ##
* The SQLite library available [here](https://www.sqlite.org/download.html).
* The Boost library available [here](http://www.boost.org).
* The fast-cpp-csv-parser by Ben Strasser, available [here](https://github.com/ben-strasser/fast-cpp-csv-parser).
* C++ compiler (such as GNU g++) with C++11 support.
* Access to the relevant folders for Default simulations (see FAQS).


## CLASS DESCRIPTIONS ##


A brief description of the important classes is given below. Some classes also contain customised exceptions for better tracing of error handling.

* The `Tree` class.
	- The most important class!
	- Contains the main setup, run and data output routines. 
	- Setup imports the data files from csv (if necessary) and creates the in-memory objects for the storing of the coalescence tree and the spatial grid of active lineages. Setup time mostly depends on the size of the csv file being imported.
	- Run continually loops over sucessive coalesence, move or speciation events until all individuals have speciated or coalesced. This is where the majority of the simulation time will be, and is mostly dependent on the number of individuals, speciation rate and size of the spatial grid.
	- At the end of the simulation, the sqlCreate() routine will generate the in-memory SQLite database for storing the coalescent tree. It can run multiple times if multiple speciation rates are required. createAndOutputData() will then be called to create a small csv file containing important information, and output the SQLite database to file if required.
* The `TreeNode` class
	- Contains a single record of a node on the phylogenetic tree, to be used in reassembling the tree structure at the end of the simulation.
	- Operations are mostly basic getters and setters, with functionality called from higher-level functions.
	- An array of treenodes makes up the `data` object in `Tree`.
* The `DataPoint` class
	- Contains a single record of the location of a lineage.
	- An array of datapoints makes up the `active` object in `Tree`.
	- `endactive` refers to the number of lineages currently being simulated. After each coalescence or speciation event this will decrease.
* The `NRrand` class
	- Contains the random number generator, as written by James Rosindell (j.rosindell@imperial.ac.uk).
* The `Landscape` class
	- Contains the routines for importing and calling values from the map objects.
	- The `getVal()` and `runDispersal()` functions can be modified to produce altered dispersal behaviour, or alterations to the structure of the map.
* The `Matrix` and `Row` classes
	- Based on code written by James Rosindell (j.rosindell@imperial.ac.uk).
	- Handles indexing of the 2D object plus importing values from a csv file.
* The `SpeciesList` class
	- Contains the species_id_list of individuals, for application in a matrix, to essentially create a 3D array.
	- Handles the positioning of individuals in space within a grid cell.
* The `ConfigParser` class
	- Contains basic functions for importing command line arguments from a config file, providing an alternative way of setting up simulations.
* The `Community` class
	 - Provides the routines for applying different speciation rates to a phylogenetic tree, to be used either immediately after simulation within necsim, or at a later time using SpeciationCounter.cpp
	 
## KNOWN BUGS ##
* Simulations run until completion, rather than aiming for a desired number of species. This is an intentional change. Functions related to this functionality remain but are deprecated.
* Only continuous rectangular fragments are properly calculated. Other shapes must be calculated by post-processing.
* 3 fragments instead of 2 will be calculated for certain adjacent rectangular patches.

## FAQS (WIP) ##
* **How do I get started?**
	- It is recommended to use the [pycoalescence](http://pycoalescence.readthedocs.io/) package which simplifies installation of necsim, setting up and running simulations. This provides a much easier way to get started with necsim.

* **Why can’t I compile the program?**
	- This could be due to a number of reasons, most likely that you haven’t compiled with access to the lsqlite3 or boost packages. Installation and compilation differs across different systems; for most UNIX systems, compiling with the linker arguments -lsqlite3 -lboost_filesystem and -lboost_system will solve problems with the compiler not finding the sqlite or boost header file.
 	- Another option could be the potential lack of access to the fast-cpp-csv-parser by Ben Strasser, available [here](https://github.com/ben-strasser/fast-cpp-csv-parser). If use\_csv has been defined at the head of the file, try without use_csv or download the csv parser and locate the folder within your working directory at compilation.
	
	
* **Every time the program runs I get error code XXX.**
	- Check the ERROR_REF.txt file for descriptions of the files. Try running in debug mode by compiling with ```-DDEBUG``` to gain more information on the problem. Check the log output in /logs. It is most likely a problem with the set up of the map data (error checking is not yet properly implemented here).
  
## CONTACTS##
Author: **Samuel Thompson**

Contact: samuelthompson14@imperial.ac.uk - thompsonsed@gmail.com

Institution: Imperial College London and National University of Singapore

Based heavily on code by **James Rosindell**

Contact: j.rosindell@imperial.ac.uk

Institution: Imperial College London
