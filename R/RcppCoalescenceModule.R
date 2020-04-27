#' rcoalescence: Efficient spatial and non-spatial neutral ecology simulator written in C++
# @useDynLib rcoalescence, .registration=TRUE

#' @importFrom Rcpp evalCpp
#' @importFrom Rcpp setRcppClass
#' @importFrom methods new
#' @import methods
#' @import RSQLite
#' @import Rcpp
#' @import dplyr
#' @name rcoalescence
#' @author Sam Thompson
#' @description A convenience wrapper for necsim, a Neutral Ecology Coalescence SIMulator for
#' spatially-explicit neutral models based in C++.
#'
#' rcoalescence allows for large-scale simulations to be performed on fully spatial files, with
#' features to suit a wide variety of use cases.
#'
#' Four classes exist to represent non-spatial/spatial and non-protracted/protracted simulations.
#' See \link[=NeutralTreeSimulation]{here} for a full description.
#'
#' @section Installation:
#' There are several dependencies required before installation of this package will be successful.
#'
#' \subsection{Dependencies}{
#'
#' \itemize{
#'     \item R package development prerequisites, which may require additional installation steps,
#'     depending on your platform
#'     \href{https://support.rstudio.com/hc/en-us/articles/200486498}{(see here)}.
#'     \item A C++ compiler that is detectable by R.
#'     \item The devtools package.
#'     \item Gdal library for C++, available \href{http://www.gdal.org/}{here}. Installed automatically for windows.
#'     \item Sqlite3 library for C++, available \href{https://www.sqlite.org/index.html}{here}.
#'
#'     \item Several R packages are dependencies, including Rcpp. These should be installed
#'     automatically.
#' }
#' }
#'
#' \subsection{Install method}{
#'
#' First, make sure devtools and the C++ library dependencies are installed. Then run,
#' \code{library(devtools)} and \code{install_github("thompsonsed/rcoalescence")}
#'
#' }
#'
#' @section Installation issues:
#'
#' \subsection{Cannot find C++ compiler}{
#'
#' For Linux systems, the R development package is required. Use
#' \code{sudo apt-get install r-base-dev} for Ubuntu-based systems. For other distributions, check
#' the documentation for installing core software development utilities required for R package
#' development.
#'
#' On MacOS, this requires that XCode is installed properly, and the compiler is detectable by R.
#' See \href{https://support.rstudio.com/hc/en-us/articles/200486498}{here} for additional info.
#'
#' Windows is currently not supported.
#' }
#'
#' \subsection{Missing gdal requirements}{
#'
#' If you get an error about gdal not being detected, check that
#' gdal is properly installed. If the error indicates missing gdal_priv.h/cpl_conv.h, consider
#' compiling gdal directly from source, as many package management systems do not properly install
#' these files. If you know the directory containing gdal, supply it manually to using PKG_LIBS=/dir/
#' during installation.
#' }
#'
#' \subsection{Missing boost requirements}{
#'
#' Boost should not be essential for package functionality. Instead, you can have a C++17 compiler,
#' or depend on the included file system methodsEnsure that boost has been properly compiled.
#' If you still require boost support, try installing the package BH from CRAN - if this
#' works then you should be fine. Alternatively, if you know the directory containing boost, supply
#' it manually to using PKG_LIBS=/dir/ during installation.
#' }
#'
#' \subsection{Still facing issues}{
#'
#' Try downloading the package from github manually and use \code{install_local()}.
#' If you have autotools installed, run \code{autoreconf} and \code{autoconf} from the before
#' attempting package
#' }
#'
#' @section Starting a simulation:
#'
#' Methods are contained within the relevant classes \link[=NeutralTreeSimulation]{(see here)} for
#' setting up, running and analysing simulations. The general process is
#' \itemize{
#'     \item Create a new \link{TreeSimulation}, \link{SpatialTreeSimulation},
#'           \link{ProtractedTreeSimulation} or \link{ProtractedSpatialTreeSimulation} object,
#'           depending on use case.
#'     \item Set the simulation parameters and add additional historical maps.
#'     \item Run the simulation (generally the most time-consuming step).
#'     \item Apply speciation rates and other requirements post simulation.
#'     \item Get the species richness from the simulation.
#'     \item Optionally write results to an SQLite database for further manual analysis.
#'     \item If written to an output database, acquire the desired biodiversity metrics, such as
#'     species locations or species abundances.
#'
#' }
#'
#' See examples in the \link[=NeutralTreeSimulation]{relevant classes}.
#'
#' @docType package
#'
"_PACKAGE"

#' Simulate neutral models on a variety of landscape conformations
#'
#' @name NeutralTreeSimulation
#' @description Four classes are provided here, \link{TreeSimulation}, \link{SpatialTreeSimulation},
#' \link{ProtractedTreeSimulation} and \link{ProtractedSpatialTreeSimulation}, for
#' non-spatial/spatial and non-protracted/protracted simulations. The process for setting up,
#' running and analysing simulations are the same, as described below.
#'
#' @section Non-spatial models:
#' Run non-spatial neutral models. Here no map files are required, simply set up the
#' simulation and run. Dispersal parameters (sigma, tau, etc) are not relevant for this class.
#'
#' @section Spatial models:
#' Run spatial neutral models. You must provide map files which provide the spatial density for the
#' simulation. Additional maps are optional.
#'
#' @section Protracted non-spatial models:
#' Run non-spatial neutral models with protracted speciation. Here no map files are required, but
#' the protracted speciation parameters should be specified.
#'
#' @section Protracted spatial models:
#' Runs spatially-explicit neutral models with protracted speciation. Requires both the protracted
#' speciation parameters and map files defining spatial density.
#'
#' @section Simulation parameters:
#' These parameters are required at simulation time. The vast majority have defaults,
#' which are defined in setSimulationParameters()
#' * *task*: the task reference number
#' * *seed*: the seed for setting random number generation
#' * *output_directory*: the path to the folder for storing simulation files
#' * *max_time*: the maximum number of seconds to simulate for before pausing
#' * *desired_specnum*: the desired number of species to aim for
#' * *times_list*: list of temporal sampling points
#' * *uses_logging*: if true, all outputs are written to console
#' * *deme*: the number of individuals per cell
#' * *deme_sample*: the global sampling proportion
#' * *min_speciation_rate*: the minimum speciation rate for the simulation
#'
#' @section Spatial parameters:
#' These include dispersal parameters and maps files
#' * *sigma*: mean dispersal: distance for a normally-distributed dispersal kernel
#' * *dispersal_method*: the method of dispersal from "normal", "fat-tailed", "uniform" and "norm-uniform"
#' * *tau*: the tau parameter for a fat-tailed dispersal kernel
#' * *m_prob*: the probability of dispersing uniformly for the norm-uniform dispersal kernel
#' * *cutoff*: the maximum distance dispersal occurs for in a uniform dispersal kernel
#' * *dispersal_relative_cost*: the relative cost of moving through non-forest
#' * *restrict_self*: if true, prevents individuals from dispersing from their own cell
#' * *landscape_type*: type of landscape from "closed", "infinite" and "tiled"
#' * *dispersal_file*: a map of dispersal probabilities
#' * *reproduction_file*: a map of reproduction probabilies
#' * *fine_map_file*: fine resolution density map
#' * *coarse_map_file*: coarse resolution density map
#' * *sample_mask_file*: spatial sampling mask
#' * *grid_x_size*: x dimension of the grid
#' * *grid_y_size*: y dimension of the grid
#' * *sample_x_size*: x dimension of the sample mask
#' * *sample_y_size*: y dimension of the sample mask
#' * *sample_x_offset*: the x offset of the sample mask from the grid
#' * *sample_y_offset*: the y offset of the sample mask from the grid
#' * *fine_map_x_size*: the x dimension of the fine map
#' * *fine_map_y_size*: the y dimension of the fine map
#' * *fine_map_x_offset*: the x offset of the fine map from the sample mask
#' * *fine_map_y_offset*: the y offset of the fine map from the sample mask
#' * *coarse_map_x_size*: the x dimension of the coarse map
#' * *coarse_map_y_size*: the y dimension of the coarse map
#' * *coarse_map_x_offset*: the x offset of the coarse map from the fine map
#' * *coarse_map_y_offset*: the y offset of the coarse map from the fine map
#' * *coarse_map_scale*: the relative scale of the coarse map
#' * *uses_spatial_sampling*: if true, the sample mask defines relative sampling proportions across the map
#' * *historical_fine_map*: the historical fine map file
#' * *historical_coarse_map*: the historical coarse map file
#' * *gen_since_historical*: the number of generations since the historical state
#' * *habitat_change_rate*: the rate of change to the historical map
#'
#' @section Protracted speciation parameters:
#' These parameters should match pairwise, i.e. the nth item in min_speciation_gens should go with
#' the nth item in max_speciation_gens
#' * *min_speciation_gens*: list of minimum number of generations required before speciation
#' * *max_speciation_gens*: list of maximum number of generations required before speciation
#'
#' @section Metacommunity parameters:
#' These parameters are used for rebuilding the coalescence tree using a metacommunity, calculated
#' using either a non-spatial simulation, an analytical approximation, or an external database
#' which contains a community to draw individuals from.
#' * *metacommunity_option*: one of "simulated", "analytical" or a path to an external database
#' * *metacommunity_size*: the number of individuals in the metacommunity
#' * *metacommunity_speciation_rate*: the effective speciation rate of the metacommunity
#' * *external_reference*: an external reference for the Metacommunity.
#'
#' @section Post-simulation parameters:
#' These are for rebuilding the coalescence tree under different conditions.
#' * *output_file*: the directory to output to, defaults to "none"
#' * *use_spatial*: if true, records full spatial locations of all individuals. Default=FALSE
#' * *sample_file*: supply a mask for defining spatial sampling. Default="null"
#' * *use_fragments*: supply a file containing fragment coordinates, or TRUE to let program calculate fragments
#' * *speciation_rates*: list of speciation rates to apply
#' * *times_list*: list of times to calculate coalescence tree at
#' * *min_speciation_gens*: list of the minimum number of generations required before speciation
#' * *max_speciation_gens*: list of the maximum number of generations required before speciation
#'
#' @md
#' @example inst/extdata/examples_nonspatial.R
#' @example inst/extdata/examples_spatial.R
#' @example inst/extdata/examples_nonspatial_protracted.R
#' @example inst/extdata/examples_spatial_protracted.R
#' @example inst/extdata/examples_metacommunity.R
NULL

Rcpp::loadModule("coalescenceModule", TRUE)

#' Non-spatial neutral simulations
#' @name TreeSimulation
#' @export TreeSimulation
#' @description Run non-spatial neutral models. Here no map files are required, simply set up the
#' simulation and run. Dispersal parameters (sigma, tau, etc) are not relevant for this class.
#'
#' @section Alternative classes:
#' For alternative classes providing similar functionality, see \link[=NeutralTreeSimulation]{here}.
#' @inheritSection NeutralTreeSimulation Simulation parameters
#' @inheritSection NeutralTreeSimulation Post-simulation parameters
#' @example inst/extdata/examples_nonspatial.R
#' @example inst/extdata/examples_metacommunity.R

TreeSimulation <- setRcppClass(
  "TreeSimulation",
  "TreeSimulation",
  module = "coalescenceModule",
  fields = list(output_database = "character"),
  methods = list(
    setInitialSimulationParameters = function(task,
                                              seed,
                                              min_speciation_rate,
                                              output_directory =
                                                "output",
                                              max_time = 3600,
                                              desired_specnum = 1,
                                              times_list = c(0.0),
                                              uses_logging = NA,
                                              deme = 1,
                                              deme_sample = 1.0) {
      "Sets the initial parameters for the simulation."
      if (!is.na(uses_logging)) {
        setLoggingMode(uses_logging)
      }
      ._setKeyParameters(
        task,
        seed,
        output_directory,
        max_time,
        desired_specnum,
        times_list
      )
      ._setDeme(deme, deme_sample)
      ._setMinSpeciationRate(min_speciation_rate)
    },

    setSpeciationParameters = function(speciation_rates,
                                       metacommunity_option = NA,
                                       metacommunity_size = NA,
                                       metacommunity_speciation_rate = NA,
                                       metacommunity_external_reference = NA) {
      "Sets the speciation parameters for applying post-simulation."
      if (is.vector(speciation_rates)) {
        for (spec_rate in speciation_rates)
        {
          ._addSpeciationRate(spec_rate)
        }
        if (length(speciation_rates) > 1) {
          ._setMultipleOutput(TRUE)
        }
      }
      else {
        ._addSpeciationRate(spec_rate)
        ._setMultipleOutput(FALSE)
      }
      if (!anyNA(metacommunity_option) |
        !anyNA(metacommunity_size) |
        !anyNA(metacommunity_speciation_rate) |
        !anyNA(metacommunity_external_reference)) {
        addMetacommunityParameters(
          metacommunity_option,
          metacommunity_size,
          metacommunity_speciation_rate,
          metacommunity_external_reference
        )
      }
    },
    addMetacommunityParameters = function(metacommunity_option =
                                            NA,
                                          metacommunity_size =
                                            NA,
                                          metacommunity_speciation_rate =
                                            NA,
                                          metacommunity_external_reference =
                                            NA) {
      if (!anyNA(metacommunity_option) | !anyNA(metacommunity_size) |
        !anyNA(metacommunity_speciation_rate) |
        !anyNA(metacommunity_external_reference)) {
        if (anyNA(metacommunity_external_reference)) {
          if (length(metacommunity_option) != length(metacommunity_size) |
            length(metacommunity_speciation_rate) != length(metacommunity_size)) {
            stop("Metacommunity parameter vectors must be of equal size.")
          }
          metacommunity_external_reference <- rep(0, length(metacommunity_option))
        }
        else {
          if (length(metacommunity_external_reference) !=
            length(metacommunity_option)) {
            stop("Metacommunity parameter vectors must be of equal size.")
          }
          if (anyNA(metacommunity_size) |
            anyNA(metacommunity_speciation_rate)) {
            metacommunity_size <- rep(0, length(metacommunity_external_reference))
            metacommunity_speciation_rate <- rep(0, length(metacommunity_external_reference))
          }
        }
        if (length(metacommunity_option) != length(metacommunity_size) |
          length(metacommunity_speciation_rate) != length(metacommunity_option) |
          length(metacommunity_external_reference) != length(metacommunity_option)) {
          stop("Metacommunity parameter vectors must be of equal size.")
        }
        if (is.vector(metacommunity_option) &
          length(metacommunity_option) > 0) {
          for (i in seq(1, length(metacommunity_option)))
          {
            ._addMetacommunityParameters(
              metacommunity_size[i],
              metacommunity_speciation_rate[i],
              metacommunity_option[i],
              metacommunity_external_reference[i]
            )
          }
          ._setMultipleOutput(TRUE)
        }
      }
      else {
        stop(
          "Must supply metacommunity parameters as either option + external reference or
            size + speciation rate."
        )
      }
    },



    applySpeciationRates = function(speciation_rates = NA,
                                    output_file = "none",
                                    times_list = c(0.0),
                                    metacommunity_option = NA,
                                    metacommunity_size = NA,
                                    metacommunity_speciation_rate = NA,
                                    metacommunity_external_reference = NA) {
      "Applies the provided speciation parameters to the simulation."
      if (!anyNA(speciation_rates)) {
        setSpeciationParameters(
          speciation_rates,
          metacommunity_option,
          metacommunity_size,
          metacommunity_speciation_rate,
          metacommunity_external_reference
        )
      }
      ._applySpeciationRates(
        output_file, FALSE, "null",
        "F", times_list
      )
    },


    setSimulationParameters = function(task,
                                       seed,
                                       min_speciation_rate,
                                       output_directory =
                                         "output",
                                       max_time = 3600,
                                       desired_specnum = 1,
                                       times_list = c(0.0),
                                       uses_logging = NA,
                                       deme = 1,
                                       deme_sample =
                                         1.0) {
      "Sets simulation parameters"
      setInitialSimulationParameters(
        task,
        seed,
        min_speciation_rate,
        output_directory,
        max_time,
        desired_specnum,
        times_list,
        uses_logging,
        deme,
        deme_sample
      )
      setup()
    },


    setOutputDatabase = function(output_db) {
      "Checks the output database exists and has been created properly,
        if it has been set."
      if (is.na(output_db) |
        output_db == "not_set") {
        stop("Output database has not been set.")
      }
      else {
        .self$output_database <- output_db
        checkOutputDatabaseExists()
      }
    },

    checkOutputDatabase = function() {
      "Checks if the output database is set"
      if (length(output_database) == 0) {
        return(FALSE)
      }
      if (is.na(output_database) |
        output_database == "none" |
        output_database == "null") {
        return(FALSE)
      }
      return(TRUE)
    },

    checkOutputDatabaseExists = function() {
      "Checks that the output database exists"
      if (!checkOutputDatabase()) {
        stop("Output database has not been set.")
      }
      if (!file.exists(output_database)) {
        stop(paste("Output database does not exist at", output_database))
      }
    },

    output = function() {
      "Outputs the biodiversity data to an sql database."
      if (!checkOutputDatabase()) {
        createDefaultOutputDatabase()
      }
      ._output()
      setOutputDatabase(._getSQLDatabase())
    },

    createDefaultOutputDatabase = function() {
      "Creates the default output database string"
      .self$output_database <-
        file.path(
          getOutputDirectory(),
          paste("data_", getJobType(), "_",
            getSeed(), ".db",
            sep =
              ""
          )
        )
    },

    getSimulationParameters = function() {
      "Returns a data frame containing all the simulation parameters."
      checkOutputDatabaseExists()
      conn <-
        dbConnect(SQLite(), output_database)
      species_locations <- dbGetQuery(
        conn,
        "SELECT seed, job_type, output_dir, speciation_rate ,sigma, tau,deme,sample_size,
          max_time,dispersal_relative_cost,min_num_species,habitat_change_rate,
          gen_since_historical,time_config_file,coarse_map_file,coarse_map_x,coarse_map_y,
          coarse_map_x_offset,coarse_map_y_offset,coarse_map_scale,fine_map_file,fine_map_x,
          fine_map_y,fine_map_x_offset,fine_map_y_offset,sample_file,grid_x,grid_y,sample_x,
          sample_y,sample_x_offset,sample_y_offset,historical_coarse_map,historical_fine_map,
          sim_complete,dispersal_method,m_probability,cutoff,restrict_self,landscape_type,
          protracted,min_speciation_gen,max_speciation_gen,dispersal_map FROM SIMULATION_PARAMETERS"
      )
      dbDisconnect(conn)
      return(species_locations)
    },

    getCommunityReferences = function() {
      "Returns a data frame containing all the community references and their
        parameter sets."
      checkOutputDatabaseExists()
      conn <-
        dbConnect(SQLite(), output_database)
      species_locations <- tryCatch(
        {
          dbGetQuery(
            conn,
            "SELECT reference, speciation_rate, time, fragments,
            metacommunity_reference, min_speciation_gen,
            max_speciation_gen FROM COMMUNITY_PARAMETERS"
          )
        },
        error = function(e) {
          return(
            dbGetQuery(
              conn,
              "SELECT reference, speciation_rate, time, fragments,
              metacommunity_reference FROM COMMUNITY_PARAMETERS"
            )
          )
        }
      )
      dbDisconnect(conn)
      return(species_locations)
    },

    getMetacommunityReferences = function() {
      "Returns a data frame containing all the metacommunity references and their
        parameter sets."
      checkOutputDatabaseExists()
      conn <-
        dbConnect(SQLite(), output_database)
      species_locations <- dbGetQuery(
        conn,
        "SELECT reference, speciation_rate,
          metacommunity_size, option,
          external_reference FROM
          METACOMMUNITY_PARAMETERS"
      )
      dbDisconnect(conn)
      return(species_locations)
    },

    getSpeciesLocations = function(community_reference = 1) {
      "Gets a data frame of species locations where the community reference
        matches the input."
      checkOutputDatabaseExists()
      conn <-
        dbConnect(SQLite(), output_database)
      community_reference_vector <- as.vector(community_reference)
      species_locations <-
        dbGetQuery(
          conn,
          paste0(
            "SELECT community_reference, species_id, x, y
              FROM SPECIES_LOCATIONS WHERE
              community_reference IN (",
            paste0("'", community_reference_vector, "'", collapse=", "), 
            ")"
          )
        )
      dbDisconnect(conn)
      if(length(community_reference_vector) == 1){
        return(species_locations %>% dplyr::select(-c("community_reference")))
      }
      return(species_locations)
    },

    getSpeciesAbundances = function(community_reference = 1) {
      "Gets a data frame of species abundances where the community reference
        matches the input"
      if (!checkOutputDatabase()) {
        return(._getSpeciesAbundances(community_reference))
      }
      checkOutputDatabaseExists()
      conn <-
        dbConnect(SQLite(), output_database)
      community_reference_vector <- as.vector(community_reference)
      species_locations <-
        dbGetQuery(
          conn,
          paste0(
            "SELECT community_reference, species_id, no_individuals FROM SPECIES_ABUNDANCES WHERE
            community_reference IN (",
            paste0("'", community_reference_vector, "'", collapse=", "), 
            ")"
          )
        )
      dbDisconnect(conn)
      if(length(community_reference_vector) == 1){
        return(species_locations %>% dplyr::select("species_id", "no_individuals"))
      }
      return(species_locations)
    },

    getSpeciesRichness = function(community_reference =
                                    NA) {
      "Gets the community reference from the output database, or from the
        internal object if no community reference is supplied (this will return the
        last calculated species richness)."
      if (!checkOutputDatabase()) {
        if (is.na(community_reference)) {
          return(._getLastSpeciesRichness())
        }
        return(._getSpeciesRichness(community_reference))
      }
      if (is.na(community_reference)) {
        community_reference <- 1
      }
      community_reference_vector <- as.vector(community_reference)
      checkOutputDatabaseExists()
      conn <-
        dbConnect(SQLite(), output_database)
      species_locations <-
        dbGetQuery(
          conn,
          paste0(
            "SELECT community_reference, COUNT(DISTINCT(species_id)) FROM SPECIES_ABUNDANCES WHERE
              no_individuals > 0 AND ",
            "community_reference IN (",
            paste0("'", community_reference_vector, "'", collapse=", "), 
            ") GROUP BY community_reference"
          )
        )
      dbDisconnect(conn)
      names(species_locations) <- c("community_reference", "species_richness")
      if(length(community_reference_vector) == 1){
        return(species_locations$species_richness[1])
      }
      return(species_locations)
    }
  )
)

#' Spatially-explicit neutral models
#' @name SpatialTreeSimulation
#' @export SpatialTreeSimulation
#' @description Run spatial neutral models. You must provide map files which provide the spatial
#' density for the simulation. Additional maps are optional.
#' @inheritSection TreeSimulation Alternative classes
#' @inheritSection NeutralTreeSimulation Spatial models
#' @inheritSection NeutralTreeSimulation Simulation parameters
#' @inheritSection NeutralTreeSimulation Spatial parameters
#' @inheritSection NeutralTreeSimulation Post-simulation parameters
#' @example inst/extdata/examples_spatial.R
SpatialTreeSimulation <- setRcppClass(
  "SpatialTreeSimulation",
  "SpatialTreeSimulation",
  module = "coalescenceModule",
  fields = list(output_database = "character"),
  contains = c("TreeSimulation"),
  methods = list(
    setDispersalParameters = function(sigma=1,
                                      dispersal_method = "normal",
                                      tau = 1.0,
                                      m_prob = 0.0,
                                      cutoff = 0,
                                      dispersal_relative_cost =
                                        1.0,
                                      restrict_self = FALSE,
                                      landscape_type = "closed",
                                      dispersal_file = "none",
                                      reproduction_file = "null") {
      "Sets the dispersal parameters for the simulation"
      ._setDispersalParameters(
        sigma,
        dispersal_method,
        tau,
        m_prob,
        cutoff,
        dispersal_relative_cost,
        restrict_self,
        landscape_type,
        dispersal_file,
        reproduction_file
      )
    },


    setMapParameters = function(fine_map_file = "null",
                                coarse_map_file = "none",
                                sample_mask_file = "null",
                                grid_x_size = NA,
                                grid_y_size = NA,
                                sample_x_size = 0,
                                sample_y_size = 0,
                                sample_x_offset = 0,
                                sample_y_offset = 0,
                                fine_map_x_size = 0,
                                fine_map_y_size = 0,
                                fine_map_x_offset = 0,
                                fine_map_y_offset = 0,
                                coarse_map_x_size = 0,
                                coarse_map_y_size = 0,
                                coarse_map_x_offset = 0,
                                coarse_map_y_offset = 0,
                                coarse_map_scale = 1,
                                deme = 1,
                                deme_sample = 1.0,
                                uses_spatial_sampling = FALSE) {
      "Sets the map parameters for the simulation"
      if (sample_mask_file == "null") {
        sample_x_size <- fine_map_x_size
        sample_y_size <- fine_map_y_size
        if (is.na(grid_x_size) | is.na(grid_y_size)) {
          grid_x_size <- fine_map_x_size
          grid_y_size <- fine_map_y_size
        }
      }
      else {
        if (is.na(grid_x_size) | is.na(grid_y_size)) {
          grid_x_size <- sample_x_size
          grid_y_size <- sample_y_size
        }
      }
      ._setMapParameters(
        fine_map_file,
        coarse_map_file,
        sample_mask_file,
        grid_x_size,
        grid_y_size,
        sample_x_size,
        sample_y_size,
        sample_x_offset,
        sample_y_offset,
        fine_map_x_size,
        fine_map_y_size,
        fine_map_x_offset,
        fine_map_y_offset,
        coarse_map_x_size,
        coarse_map_y_size,
        coarse_map_x_offset,
        coarse_map_y_offset,
        coarse_map_scale,
        deme,
        deme_sample,
        uses_spatial_sampling
      )
    },


    setHistoricalMapParameters = function(historical_fine_map =
                                            "none",
                                          historical_coarse_map =
                                            "none",
                                          gen_since_historical =
                                            100000000,
                                          habitat_change_rate =
                                            0.0) {
      "Sets the historical map parameters for the simulation."
      ._setHistoricalMapParameters(
        historical_fine_map,
        historical_coarse_map,
        gen_since_historical,
        habitat_change_rate
      )
    },

    addHistoricalMap = function(historical_fine_map,
                                historical_coarse_map = "none",
                                gen_since_historical = 1,
                                habitat_change_rate = 0.0) {
      "Adds a historical map to the list of historical maps to use."
      ._addHistoricalMap(
        historical_fine_map,
        historical_coarse_map,
        gen_since_historical,
        habitat_change_rate
      )
    },

    applySpeciationRates = function(speciation_rates = NA,
                                    output_file = "none",
                                    use_spatial = FALSE,
                                    sample_file = "null",
                                    use_fragments = FALSE,
                                    times_list = c(0.0),
                                    metacommunity_option = NA,
                                    metacommunity_size = NA,
                                    metacommunity_speciation_rate =
                                      NA,
                                    metacommunity_external_reference =
                                      NA) {
      "Applies the provided speciation parameters to the simulation"
      if (is.logical(use_fragments)) {
        use_fragments <- substr(as.character(use_fragments), 1, 1)
      }
      if (!anyNA(speciation_rates)) {
        setSpeciationParameters(
          speciation_rates,
          metacommunity_option,
          metacommunity_size,
          metacommunity_speciation_rate,
          metacommunity_external_reference
        )
      }
      ._applySpeciationRates(
        output_file,
        use_spatial,
        sample_file,
        use_fragments,
        times_list
      )
    },


    setSimulationParameters = function(task,
                                       seed,
                                       min_speciation_rate,
                                       sigma=1,
                                       output_directory = "output",
                                       max_time = 3600,
                                       desired_specnum = 1,
                                       times_list = c(0.0),
                                       uses_logging = NA,
                                       dispersal_method = "normal",
                                       tau = 1.0,
                                       m_prob = 0.0,
                                       cutoff = 0,
                                       dispersal_relative_cost =
                                         1.0,
                                       restrict_self = FALSE,
                                       landscape_type = "closed",
                                       dispersal_file = "none",
                                       reproduction_file = "none",
                                       fine_map_file = "null",
                                       coarse_map_file = "none",
                                       sample_mask_file = "null",
                                       grid_x_size = NA,
                                       grid_y_size = NA,
                                       sample_x_size = 0,
                                       sample_y_size = 0,
                                       sample_x_offset = 0,
                                       sample_y_offset = 0,
                                       fine_map_x_size = 0,
                                       fine_map_y_size = 0,
                                       fine_map_x_offset = 0,
                                       fine_map_y_offset = 0,
                                       coarse_map_x_size = 0,
                                       coarse_map_y_size = 0,
                                       coarse_map_x_offset =
                                         0,
                                       coarse_map_y_offset =
                                         0,
                                       coarse_map_scale = 1,
                                       deme = 1,
                                       deme_sample = 1.0,
                                       uses_spatial_sampling = FALSE,
                                       historical_fine_map =
                                         "none",
                                       historical_coarse_map =
                                         "none",
                                       gen_since_historical =
                                         100000000,
                                       habitat_change_rate =
                                         0.0) {
      "Sets all simulation parameters"
      setInitialSimulationParameters(
        task,
        seed,
        min_speciation_rate,
        output_directory,
        max_time,
        desired_specnum,
        times_list,
        uses_logging
      )
      setDispersalParameters(
        sigma,
        dispersal_method,
        tau,
        m_prob,
        cutoff,
        dispersal_relative_cost,
        restrict_self,
        landscape_type,
        dispersal_file,
        reproduction_file
      )
      setMapParameters(
        fine_map_file,
        coarse_map_file,
        sample_mask_file,
        grid_x_size,
        grid_y_size,
        sample_x_size,
        sample_y_size,
        sample_x_offset,
        sample_y_offset,
        fine_map_x_size,
        fine_map_y_size,
        fine_map_x_offset,
        fine_map_y_offset,
        coarse_map_x_size,
        coarse_map_y_size,
        coarse_map_x_offset,
        coarse_map_y_offset,
        coarse_map_scale,
        deme,
        deme_sample,
        uses_spatial_sampling
      )
      setHistoricalMapParameters(
        historical_fine_map,
        historical_coarse_map,
        gen_since_historical,
        habitat_change_rate
      )
      setup()
    }
  )
)

#' Protracted non-spatial neutral models
#' @name ProtractedTreeSimulation
#' @export ProtractedTreeSimulation
#' @description Run non-spatial neutral models with protracted speciation. Here no map files are
#' required, but the protracted speciation parameters should be specified.
#' @inheritSection TreeSimulation Alternative classes
#' @inheritSection NeutralTreeSimulation Protracted non-spatial models
#' @inheritSection NeutralTreeSimulation Simulation parameters
#' @inheritSection NeutralTreeSimulation Protracted speciation parameters
#' @inheritSection NeutralTreeSimulation Post-simulation parameters
#' @example inst/extdata/examples_nonspatial_protracted.R
ProtractedTreeSimulation <- setRcppClass(
  "ProtractedTreeSimulation",
  "ProtractedTreeSimulation",
  module = "coalescenceModule",
  fields = list(output_database = "character"),
  contains = c("TreeSimulation"),
  methods = list(
    setSimulationParameters = function(task,
                                       seed,
                                       min_speciation_rate,
                                       output_directory =
                                         "output",
                                       max_time = 3600,
                                       desired_specnum = 1,
                                       times_list = c(0.0),
                                       uses_logging = NA,
                                       deme = 1,
                                       deme_sample =
                                         1.0,
                                       min_speciation_gen =
                                         NA,
                                       max_speciation_gen =
                                         NA) {
      "Sets all protracted simulation parameters"
      setInitialSimulationParameters(
        task,
        seed,
        min_speciation_rate,
        output_directory,
        max_time,
        desired_specnum,
        times_list,
        uses_logging,
        deme,
        deme_sample
      )
      if (!is.na(min_speciation_gen) &
        !is.na(max_speciation_gen)) {
        ._setSimulationProtractedParameters(
          TRUE, min_speciation_gen,
          max_speciation_gen
        )
      }
      setup()
    },

    applySpeciationRates = function(speciation_rates,
                                    output_file = "none",
                                    times_list = c(0.0),
                                    min_speciation_gens =
                                      c(0.0),
                                    max_speciation_gens =
                                      c(0.0),
                                    metacommunity_option =
                                      NA,
                                    metacommunity_size =
                                      NA,
                                    metacommunity_speciation_rate =
                                      NA,
                                    metacommunity_external_reference =
                                      NA) {
      "Applies the provided speciation parameters"
      if (!anyNA(speciation_rates)) {
        setSpeciationParameters(
          speciation_rates,
          metacommunity_option,
          metacommunity_size,
          metacommunity_speciation_rate,
          metacommunity_external_reference
        )
      }
      if (!anyNA(min_speciation_gens) &
        !anyNA(max_speciation_gens)) {
        if (is.vector(min_speciation_gens) & is.vector(max_speciation_gens)) {
          if (length(min_speciation_gens) != length(max_speciation_gens)) {
            stop("Lengths of protracted parameter vectors must be equal.")
          }
          for (i in seq(1, length(min_speciation_gens)))
          {
            ._addProtractedParameters(
              min_speciation_gens[i],
              max_speciation_gens[i]
            )
          }
          ._setMultipleOutput(TRUE)
        }
        else {
          ._addProtractedParameters(min_speciation_gens, max_speciation_gens)
        }
      }
      ._applySpeciationRates(output_file, FALSE, "null", "F", times_list)
    },


    setSimulationParameters = function(task,
                                       seed,
                                       speciation_rate,
                                       output_directory =
                                         "output",
                                       max_time = 3600,
                                       desired_specnum = 1,
                                       times_list = c(0.0),
                                       uses_logging = NA,
                                       min_speciation_gen = 0.0,
                                       max_speciation_gen = 0.0,
                                       deme = 1,
                                       deme_sample =
                                         1.0) {
      "Sets all simulation parameters"
      setInitialSimulationParameters(
        task,
        seed,
        speciation_rate,
        output_directory,
        max_time,
        desired_specnum,
        times_list,
        uses_logging,
        deme,
        deme_sample
      )
      ._setSimulationProtractedParameters(
        TRUE, min_speciation_gen,
        max_speciation_gen
      )
      setup()
    }
  )
)

#' Protracted spatially-explicit neutral models
#' @name ProtractedSpatialTreeSimulation
#' @export ProtractedSpatialTreeSimulation
#' @description  Runs spatially-explicit neutral models with protracted speciation. Requires both
#' the protracted speciation parameters and map files defining spatial density.
#' @inheritSection TreeSimulation Alternative classes
#' @inheritSection NeutralTreeSimulation Protracted spatial models
#' @inheritSection NeutralTreeSimulation Simulation parameters
#' @inheritSection NeutralTreeSimulation Spatial parameters
#' @inheritSection NeutralTreeSimulation Protracted speciation parameters
#' @inheritSection NeutralTreeSimulation Post-simulation parameters
#' @example inst/extdata/examples_spatial_protracted.R
ProtractedSpatialTreeSimulation <-setRcppClass(
  "ProtractedSpatialTreeSimulation",
  "ProtractedSpatialTreeSimulation",
  module = "coalescenceModule",
  fields = list(output_database = "character"),
  contains = c(
    "ProtractedTreeSimulation",
    "SpatialTreeSimulation",
    "TreeSimulation"
  ),
  methods = list(
    setSimulationParameters = function(task,
                                       seed,
                                       min_speciation_rate,
                                       sigma=1,
                                       output_directory = "output",
                                       max_time = 3600,
                                       desired_specnum = 1,
                                       times_list = c(0.0),
                                       uses_logging = NA,
                                       dispersal_method = "normal",
                                       tau = 1.0,
                                       m_prob = 0.0,
                                       cutoff = 0,
                                       dispersal_relative_cost = 1.0,
                                       restrict_self = FALSE,
                                       landscape_type = "closed",
                                       dispersal_file = "none",
                                       reproduction_file = "none",
                                       fine_map_file = "null",
                                       coarse_map_file = "none",
                                       sample_mask_file = "null",
                                       grid_x_size = NA,
                                       grid_y_size = NA,
                                       sample_x_size = 0,
                                       sample_y_size = 0,
                                       sample_x_offset = 0,
                                       sample_y_offset = 0,
                                       fine_map_x_size = 0,
                                       fine_map_y_size = 0,
                                       fine_map_x_offset = 0,
                                       fine_map_y_offset = 0,
                                       coarse_map_x_size = 0,
                                       coarse_map_y_size = 0,
                                       coarse_map_x_offset = 0,
                                       coarse_map_y_offset = 0,
                                       coarse_map_scale = 1,
                                       deme = 1,
                                       deme_sample = 1.0,
                                       uses_spatial_sampling = FALSE,
                                       historical_fine_map = "none",
                                       historical_coarse_map = "none",
                                       gen_since_historical = 100000000,
                                       habitat_change_rate = 0.0,
                                       min_speciation_gen = 0.0,
                                       max_speciation_gen = 0.0) {
      "Sets all simulation parameters"
      setInitialSimulationParameters(
        task,
        seed,
        min_speciation_rate,
        output_directory,
        max_time,
        desired_specnum,
        times_list,
        uses_logging
      )
      ._setSimulationProtractedParameters(
        TRUE, min_speciation_gen,
        max_speciation_gen
      )
      setDispersalParameters(
        sigma,
        dispersal_method,
        tau,
        m_prob,
        cutoff,
        dispersal_relative_cost,
        restrict_self,
        landscape_type,
        dispersal_file,
        reproduction_file
      )
      setMapParameters(
        fine_map_file,
        coarse_map_file,
        sample_mask_file,
        grid_x_size,
        grid_y_size,
        sample_x_size,
        sample_y_size,
        sample_x_offset,
        sample_y_offset,
        fine_map_x_size,
        fine_map_y_size,
        fine_map_x_offset,
        fine_map_y_offset,
        coarse_map_x_size,
        coarse_map_y_size,
        coarse_map_x_offset,
        coarse_map_y_offset,
        coarse_map_scale,
        deme,
        deme_sample,
        uses_spatial_sampling
      )
      setHistoricalMapParameters(
        historical_fine_map,
        historical_coarse_map,
        gen_since_historical,
        habitat_change_rate
      )
      setup()
    },

    applySpeciationRates = function(speciation_rates,
                                    output_file = "none",
                                    use_spatial = FALSE,
                                    sample_file = "null",
                                    use_fragments = FALSE,
                                    times_list = c(0.0),
                                    min_speciation_gens = c(0.0),
                                    max_speciation_gens = c(0.0),
                                    metacommunity_option = NA,
                                    metacommunity_size = NA,
                                    metacommunity_speciation_rate = NA,
                                    metacommunity_external_reference = NA) {
      "Applies the provided speciation parameters to the simulation"
      if (is.logical(use_fragments)) {
        use_fragments <- substr(as.character(use_fragments), 1, 1)
      }
      if (!anyNA(speciation_rates)) {
        setSpeciationParameters(
          speciation_rates,
          metacommunity_option,
          metacommunity_size,
          metacommunity_speciation_rate,
          metacommunity_external_reference
        )
      }
      if (!anyNA(min_speciation_gens) &
        !anyNA(max_speciation_gens)) {
        if (is.vector(min_speciation_gens) & is.vector(max_speciation_gens)) {
          if (length(min_speciation_gens) != length(max_speciation_gens)) {
            stop("Lengths of protracted parameter vectors must be equal.")
          }
          for (i in seq(1, length(min_speciation_gens)))
          {
            ._addProtractedParameters(
              min_speciation_gens[i],
              max_speciation_gens[i]
            )
          }
          ._setMultipleOutput(TRUE)
        }
        else {
          ._addProtractedParameters(min_speciation_gens, max_speciation_gens)
        }
      }
      ._applySpeciationRates(
        output_file,
        use_spatial,
        sample_file,
        use_fragments,
        times_list
      )
    }
  )
)
