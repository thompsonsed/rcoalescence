# Tests that coalescence simulations run as expected

context("Basic spatial coalescence simulations")
test_that("Basic simulation on a null landscape completes", {
  tmp <- new("SpatialTreeSimulation")
  tmp$setInitialSimulationParameters(
    task = 0,
    seed =  0,
    min_speciation_rate = 0.1,
    output_directory = "default",
    max_time =  10,
    desired_specnum = 1,
    times_list = c(0.0)
  )
  tmp$setDispersalParameters(sigma = 2)
  tmp$setHistoricalMapParameters()
  tmp$setMapParameters(fine_map_x_size = 10, fine_map_y_size = 10)
  tmp$setup()
  expect_equal(TRUE, tmp$runSimulation())
})

test_that("Basic simulation produces the expected number of individuals in edge cases",
          {
            tmp <- SpatialTreeSimulation$new()
            tmp$setSimulationParameters(
              task = 12,
              seed = 1,
              min_speciation_rate = 0.000000000001,
              deme = 1,
              sigma = 2,
              fine_map_x_size = 10,
              fine_map_y_size = 10,
              uses_logging = FALSE
            )
            expect_equal(TRUE, tmp$runSimulation())
            tmp$applySpeciationRates(speciation_rates = 0.000000000001)
            expect_equal(1, tmp$getSpeciesRichness())
            tmp$applySpeciationRates(speciation_rates = 0.99999)
            expect_equal(100, tmp$getSpeciesRichness())
          })

test_that("Running with a basic map file works as expected", {
  tmp <- SpatialTreeSimulation$new()
  tmp$setSimulationParameters(
    seed = 9,
    task = 1,
    output_directory = "output",
    min_speciation_rate = 0.5,
    sigma = 2 * (2 ** 0.5),
    deme = 1,
    deme_sample = 0.1,
    fine_map_file = "../../inst/extdata/sample/example_fine.tif",
    fine_map_x_size = 13,
    fine_map_y_size = 13,
    coarse_map_file = "../../inst/extdata/sample/example_coarse.tif",
    coarse_map_x_size = 35,
    coarse_map_y_size = 41,
    coarse_map_x_offset = 11,
    coarse_map_y_offset = 14,
    sample_mask_file = "../../inst/extdata/sample/example_mask.tif",
    sample_x_size = 13,
    sample_y_size = 13,
    uses_logging = FALSE
  )
  expect_equal(TRUE, tmp$runSimulation())
  tmp$applySpeciationRates(speciation_rates = 0.5)
  expect_equal(1164, tmp$getSpeciesRichness())
})

test_that("Simulations correctly detect map sizes", {
  tmp <- SpatialTreeSimulation$new()
  tmp$setSimulationParameters(
    seed = 9,
    task = 1,
    output_directory = "output",
    min_speciation_rate = 0.5,
    sigma = 2 * (2 ** 0.5),
    deme = 1,
    deme_sample = 0.1,
    fine_map_file = "../../inst/extdata/sample/example_fine.tif",
    coarse_map_file = "../../inst/extdata/sample/example_coarse.tif",
    sample_mask_file = "../../inst/extdata/sample/example_mask.tif",
    uses_logging = FALSE
  )
  expect_equal(TRUE, tmp$runSimulation())
  tmp$applySpeciationRates(speciation_rates = c(0.5))
  expect_equal(1164, tmp$getSpeciesRichness())
  expect_error(tmp$getSimulationParameters())
  tmp$output()
  sim_params <- tmp$getSimulationParameters()
  expected_params <-
    structure(
      list(
        seed = 9L,
        job_type = 1L,
        output_dir = "output",
        speciation_rate = 0.5,
        sigma = 2.828427,
        tau = 1,
        deme = 1L,
        sample_size = 0.1,
        max_time = 3600L,
        dispersal_relative_cost = 1,
        min_num_species = 1L,
        habitat_change_rate = 0,
        gen_since_historical = 1e+08,
        time_config_file = "set",
        coarse_map_file = "../../inst/extdata/sample/example_coarse.tif",
        coarse_map_x = 35L,
        coarse_map_y = 41L,
        coarse_map_x_offset = 11L,
        coarse_map_y_offset = 14L,
        coarse_map_scale = 1,
        fine_map_file = "../../inst/extdata/sample/example_fine.tif",
        fine_map_x = 13L,
        fine_map_y = 13L,
        fine_map_x_offset = 0L,
        fine_map_y_offset = 0L,
        sample_file = "../../inst/extdata/sample/example_mask.tif",
        grid_x = 13L,
        grid_y = 13L,
        sample_x = 13L,
        sample_y = 13L,
        sample_x_offset = 0L,
        sample_y_offset = 0L,
        historical_coarse_map = "none",
        historical_fine_map = "none",
        sim_complete = 1L,
        dispersal_method = "normal",
        m_probability = 0,
        cutoff = 0,
        restrict_self = 0L,
        landscape_type = "closed",
        protracted = 0L,
        min_speciation_gen = 0,
        max_speciation_gen = 0,
        dispersal_map = "none"
      ),
      .Names = c(
        "seed",
        "job_type",
        "output_dir",
        "speciation_rate",
        "sigma",
        "tau",
        "deme",
        "sample_size",
        "max_time",
        "dispersal_relative_cost",
        "min_num_species",
        "habitat_change_rate",
        "gen_since_historical",
        "time_config_file",
        "coarse_map_file",
        "coarse_map_x",
        "coarse_map_y",
        "coarse_map_x_offset",
        "coarse_map_y_offset",
        "coarse_map_scale",
        "fine_map_file",
        "fine_map_x",
        "fine_map_y",
        "fine_map_x_offset",
        "fine_map_y_offset",
        "sample_file",
        "grid_x",
        "grid_y",
        "sample_x",
        "sample_y",
        "sample_x_offset",
        "sample_y_offset",
        "historical_coarse_map",
        "historical_fine_map",
        "sim_complete",
        "dispersal_method",
        "m_probability",
        "cutoff",
        "restrict_self",
        "landscape_type",
        "protracted",
        "min_speciation_gen",
        "max_speciation_gen",
        "dispersal_map"
      ),
      class = "data.frame",
      row.names = c(NA, -1L)
    )
  expect_equal(TRUE, all.equal(sim_params, expected_params))
  
})

context("Checking biodiversity metrics simulations")
test_that("Biodiversity metrics correctly stored in output database", {
  output_file = file.path("output", "data_1_10.db")
  if (file.exists(output_file))
  {
    file.remove(output_file)
  }
  tmp <- SpatialTreeSimulation$new()
  tmp$setSimulationParameters(
    seed = 10,
    task = 1,
    output_directory = "output",
    min_speciation_rate = 0.5,
    sigma = 2 * (2 ** 0.5),
    deme = 1,
    deme_sample = 0.1,
    fine_map_file = "../../inst/extdata/sample/example_fine.tif",
    coarse_map_file = "../../inst/extdata/sample/example_coarse.tif",
    sample_mask_file = "../../inst/extdata/sample/example_mask.tif",
    uses_logging = FALSE
  )
  expect_equal(TRUE, tmp$runSimulation())
  tmp$applySpeciationRates(speciation_rates = c(0.5, 0.7),
                           use_spatial = TRUE)
  tmp$output()
  # Make sure an output exists
  expect_equal(TRUE, file.exists(output_file))
  # Check community references
  community_references <-
    data.frame(matrix(
      c(1, 0.5, 0.0, 0, 0, 2, 0.7, 0.0, 0, 0),
      nrow = 2,
      ncol = 5,
      byrow = TRUE
    ))
  names(community_references) <-
    c("reference",
      "speciation_rate",
      "time",
      "fragments",
      "metacommunity_reference")
  expect_equal(TRUE,
               all.equal(tmp$getCommunityReferences(), community_references))
  # Check species abundances
  expected_abundances <-
    data.frame(matrix(c(0, 0, 1, 1, 2, 1, 3, 1, 4, 1), ncol = 2, byrow = TRUE))
  names(expected_abundances) <- c("species_id", "no_individuals")
  actual_abundances <- head(tmp$getSpeciesAbundances(1), 5)
  expect_equal(TRUE, all.equal(expected_abundances, actual_abundances))
  expected_abundances <-
    data.frame(matrix(
      c(1162, 2, 1163, 2, 1164, 2, 1165, 1, 1166, 1),
      ncol = 2,
      byrow = TRUE
    ))
  names(expected_abundances) <- c("species_id", "no_individuals")
  actual_abundances <-
    data.frame(tail(tmp$getSpeciesAbundances(2), 5), row.names = NULL)
  expect_equal(TRUE, isTRUE(all.equal(
    expected_abundances, actual_abundances
  )))
  # Check species locations
  expected_locations <-
    data.frame(matrix(
      c(1, 0, 10, 2, 0, 10, 3, 0, 10, 4, 0, 10, 5, 0, 10),
      ncol = 3,
      byrow = TRUE
    ))
  names(expected_locations) <- c("species_id", "x", "y")
  actual_locations <- head(tmp$getSpeciesLocations(1), 5)
  expect_equal(TRUE, all.equal(expected_locations, actual_locations))
  expected_locations <-
    data.frame(matrix(
      c(1150, 12, 12, 1151, 12, 12, 1152, 12, 12, 1153, 12, 12, 1154, 12, 12),
      ncol = 3,
      byrow = TRUE
    ))
  names(expected_locations) <- c("species_id", "x", "y")
  actual_locations <-
    data.frame(tail(tmp$getSpeciesLocations(2), 5), row.names = NULL)
  expect_equal(TRUE, isTRUE(all.equal(expected_locations, actual_locations)))
  # Check species richness
  expect_equal(1160, tmp$getSpeciesRichness(1))
  expect_equal(1166, tmp$getSpeciesRichness(2))
  
  
  
})
context("More complex spatial coalescence simulations")
test_that("Simulation with a single historical maps works as intended.", {
  tmp <- SpatialTreeSimulation$new()
  tmp$setSimulationParameters(
    seed = 10,
    task = 1,
    output_directory = "output",
    min_speciation_rate = 0.1,
    sigma = 2 * (2 ** 0.5),
    deme = 1,
    deme_sample = 0.1,
    fine_map_file = "../../inst/extdata/sample/example_fine.tif",
    coarse_map_file = "../../inst/extdata/sample/example_coarse.tif",
    sample_mask_file = "../../inst/extdata/sample/example_mask.tif",
    uses_logging = FALSE,
    landscape_type = "closed"
  )
  tmp$addHistoricalMap(
    historical_fine_map = "../../inst/extdata/sample/example_historical_fine.tif",
    historical_coarse_map = "../../inst/extdata/sample/example_coarse.tif",
    gen_since_historical = 1,
    habitat_change_rate = 0.5
  )
  # tmp$addHistoricalMap(historical_fine_map="../../inst/extdata/sample/example_fine.tif",
  #                      historical_coarse_map="../../inst/extdata/sample/example_coarse.tif",
  #                      gen_since_historical=10, habitat_change_rate=0.5)
  # tmp$addHistoricalMap(historical_fine_map="../../inst/extdata/sample/example_fine.tif",
  # historical_coarse_map="null", gen_since_historical=4,
  # habitat_change_rate=1.0)
  expect_equal(TRUE, tmp$runSimulation())
  tmp$applySpeciationRates(speciation_rates = c(0.1))
  expect_equal(1118, tmp$getSpeciesRichness())
})
#
test_that("Simulation with multiple historical maps works as intended.", {
  tmp <- SpatialTreeSimulation$new()
  tmp$setSimulationParameters(
    seed = 10,
    task = 1,
    output_directory = "output",
    min_speciation_rate = 0.5,
    sigma = 2 * (2 ** 0.5),
    deme = 1,
    deme_sample = 0.1,
    fine_map_file = "../../inst/extdata/sample/example_fine.tif",
    coarse_map_file = "../../inst/extdata/sample/example_coarse.tif",
    sample_mask_file = "../../inst/extdata/sample/example_mask.tif",
    uses_logging = FALSE
  )
  tmp$addHistoricalMap(
    historical_fine_map = "../../inst/extdata/sample/example_historical_fine.tif",
    historical_coarse_map = "../../inst/extdata/sample/example_coarse.tif",
    gen_since_historical = 1,
    habitat_change_rate = 0.5
  )
  tmp$addHistoricalMap(
    historical_fine_map = "null",
    historical_coarse_map = "null",
    gen_since_historical = 10,
    habitat_change_rate = 1.0
  )
  tmp$addHistoricalMap(
    historical_fine_map = "../../inst/extdata/sample/example_fine.tif",
    historical_coarse_map = "../../inst/extdata/sample/example_coarse.tif",
    gen_since_historical = 20,
    habitat_change_rate = 1.0
  )
  expect_equal(TRUE, tmp$runSimulation())
  tmp$applySpeciationRates(speciation_rates = c(0.5))
  expect_equal(1159, tmp$getSpeciesRichness())
})
#
test_that("Simulation with multiple sampling times.", {
  tmp <- SpatialTreeSimulation$new()
  tmp$setSimulationParameters(
    seed = 11,
    task = 1,
    output_directory = "output",
    min_speciation_rate = 0.5,
    sigma = 2 * (2 ** 0.5),
    deme = 1,
    deme_sample = 0.1,
    fine_map_file = "../../inst/extdata/sample/example_fine.tif",
    coarse_map_file = "../../inst/extdata/sample/example_coarse.tif",
    sample_mask_file = "../../inst/extdata/sample/example_mask.tif",
    uses_logging = FALSE,
    times_list = c(0.0, 1.0, 10.0, 20.0)
  )
  expect_equal(TRUE, tmp$runSimulation())
  tmp$applySpeciationRates(speciation_rates = c(0.5),
                           times_list = c(0.0, 1.0, 10.0, 20.0))
  tmp$output()
  community_references <-
    data.frame(matrix(
      c(
        1,
        0.5,
        0.0,
        0,
        0,
        2,
        0.5,
        1.0,
        0,
        0,
        3,
        0.5,
        10.0,
        0,
        0,
        4,
        0.5,
        20.0,
        0,
        0
      ),
      ncol = 5,
      byrow = TRUE
    ))
  names(community_references) <-
    c("reference",
      "speciation_rate",
      "time",
      "fragments",
      "metacommunity_reference")
  expect_equal(TRUE,
               all.equal(community_references, tmp$getCommunityReferences()))
  expect_equal(1159, tmp$getSpeciesRichness(1))
  expect_equal(1160, tmp$getSpeciesRichness(2))
  expect_equal(1154, tmp$getSpeciesRichness(3))
  expect_equal(1163, tmp$getSpeciesRichness(4))
})


context("Basic non-spatial coalescence simulations")
test_that("Basic simulation with deme of 100 completes", {
  tmp <- TreeSimulation$new()
  tmp$setSimulationParameters(
    task = 101,
    seed =  1,
    output_directory =  "default",
    max_time = 10,
    deme = 100,
    uses_logging = FALSE,
    min_speciation_rate = 0.1
  )
  expect_equal(TRUE, tmp$runSimulation())
  tmp$applySpeciationRates(speciation_rates = 0.1)
  expect_equal(24, tmp$getSpeciesRichness())
})

test_that("Basic simulation with deme of 100 completes and can output without application",
          {
            tmp <- TreeSimulation$new()
            tmp$setSimulationParameters(
              task = 101,
              seed =  1,
              output_directory =  "default",
              max_time = 10,
              deme = 100,
              uses_logging = FALSE,
              min_speciation_rate = 0.1
            )
            expect_equal(TRUE, tmp$runSimulation())
            tmp$output()
            expect_equal(TRUE, file.exists(tmp$output_database))
          })

test_that("Basic simulation produces 100 species with maximum speciation rate",
          {
            tmp <- TreeSimulation$new()
            tmp$setSimulationParameters(
              task = 102,
              seed = 1,
              min_speciation_rate = 0.99999,
              deme = 100,
              uses_logging = FALSE
            )
            expect_equal(TRUE, tmp$runSimulation())
            tmp$applySpeciationRates(speciation_rates = 0.999999)
            expect_equal(100, tmp$getSpeciesRichness())
          })

test_that("Basic simulation produces 1 species with minimum speciation rate",
          {
            tmp <- TreeSimulation$new()
            tmp$setSimulationParameters(
              task = 103,
              seed = 2,
              min_speciation_rate = 0.00000000001,
              deme = 100,
              uses_logging = FALSE
            )
            expect_equal(TRUE, tmp$runSimulation())
            tmp$applySpeciationRates(speciation_rates = 0.00000000001)
            expect_equal(1, tmp$getSpeciesRichness())
          })



context("Basic non-spatial protracted coalescence simulations")
test_that("Basic simulation with deme of 100 completes", {
  tmp <- ProtractedTreeSimulation$new()
  tmp$setSimulationParameters(
    task = 201,
    seed = 1,
    output_directory = "output",
    max_time = 200,
    desired_specnum = 1,
    times_list = c(0.0),
    deme = 100,
    min_speciation_rate = 0.01,
    min_speciation_gen = 2.0,
    max_speciation_gen = 4000.0,
    uses_logging = FALSE
  )
  expect_equal(TRUE, tmp$runSimulation())
  tmp$applySpeciationRates(
    speciation_rates = c(0.02, 0.03),
    min_speciation_gens = c(2.0),
    max_speciation_gens = c(4000.0)
  )
  tmp$output()
  expect_equal(6, tmp$getSpeciesRichness())
  expect_equal(6, tmp$getSpeciesRichness(1))
  expect_equal(12, tmp$getSpeciesRichness(2))
})

test_that("Basic simulation produces 1 species with high speciation generation",
          {
            tmp <- ProtractedTreeSimulation$new()
            tmp$setSimulationParameters(
              task = 202,
              seed = 3,
              min_speciation_rate = 0.9999,
              deme = 100,
              uses_logging = FALSE,
              min_speciation_gen = 1000000.0,
              max_speciation_gen = 11000000.0
            )
            expect_equal(TRUE, tmp$runSimulation())
            tmp$applySpeciationRates(
              speciation_rates = 0.9999,
              min_speciation_gens = c(1000000.0),
              max_speciation_gens = c(11000000.0)
            )
            expect_equal(1, tmp$getSpeciesRichness())
          })

test_that("Basic simulation produces 100 species with low speciation generation",
          {
            tmp <- ProtractedTreeSimulation$new()
            tmp$setSimulationParameters(
              task = 203,
              seed =  3,
              min_speciation_rate = 0.000001,
              deme = 100,
              uses_logging = FALSE,
              min_speciation_gen = 0.0,
              max_speciation_gen = 0.01
            )
            expect_equal(TRUE, tmp$runSimulation())
            tmp$applySpeciationRates(
              speciation_rates = 0.000001,
              min_speciation_gens = 0.0,
              max_speciation_gens = 0.01
            )
            expect_equal(100, tmp$getSpeciesRichness())
          })

test_that("Basic simulation produces 100 species with maximum speciation rate",
          {
            tmp <- ProtractedTreeSimulation$new()
            tmp$setSimulationParameters(
              task = 204,
              seed = 1,
              min_speciation_rate = 0.1,
              deme = 100,
              uses_logging = FALSE,
              min_speciation_gen = 1.0,
              max_speciation_gen = 30.0
            )
            expect_equal(TRUE, tmp$runSimulation())
            tmp2 <- ProtractedTreeSimulation$new()
            tmp2$setSimulationParameters(
              task = 205,
              seed = 1,
              min_speciation_rate = 0.1,
              deme = 100,
              uses_logging = FALSE,
              min_speciation_gen = 1.0,
              max_speciation_gen = 10.0
            )
            expect_equal(TRUE, tmp2$runSimulation())
            tmp3 <- ProtractedTreeSimulation$new()
            tmp3$setSimulationParameters(
              task = 206,
              seed = 1,
              min_speciation_rate = 0.1,
              deme = 100,
              uses_logging = FALSE,
              min_speciation_gen = 5.0,
              max_speciation_gen = 10.0
            )
            expect_equal(TRUE, tmp3$runSimulation())
            tmp$applySpeciationRates(
              speciation_rates = 0.1,
              min_speciation_gens = 1.0,
              max_speciation_gens = 30.0
            )
            r1 <- tmp$getSpeciesRichness()
            tmp2$applySpeciationRates(
              speciation_rates = 0.1,
              min_speciation_gens = 1.0,
              max_speciation_gens = 10.0
            )
            r2 <- tmp2$getSpeciesRichness()
            tmp3$applySpeciationRates(
              speciation_rates = 0.1,
              min_speciation_gens = 5.0,
              max_speciation_gens = 10.0
            )
            r3 <- tmp3$getSpeciesRichness()
            expect_equal(TRUE, r2 > r1)
            expect_equal(TRUE, r2 > r3)
            expect_equal(28, r2)
            expect_equal(25, r1)
            expect_equal(20, r3)
          })


test_that("Can apply multiple protracted speciation parameters", {
  tmp <- ProtractedTreeSimulation$new()
  tmp$setSimulationParameters(
    task = 204,
    seed = 1,
    min_speciation_rate = 0.1,
    deme = 100,
    uses_logging = FALSE,
    min_speciation_gen = 10.0,
    max_speciation_gen = 30.0
  )
  expect_equal(TRUE, tmp$runSimulation())
  tmp$applySpeciationRates(
    speciation_rates = c(0.1, 0.2, 0.3),
    min_speciation_gens = c(1.0, 5.0, 10.0),
    max_speciation_gens = c(10.0, 10.0, 10.0)
  )
  tmp$applySpeciationRates(
    speciation_rates = c(0.1, 0.2, 0.3),
    min_speciation_gens = c(1.0, 5.0, 10.0),
    max_speciation_gens = c(20.0, 20.0, 20.0)
  )
  tmp$output()
  crefs <- tmp$getCommunityReferences()
  sr <- c(0.1,
          0.2,
          0.3,
          0.1,
          0.2,
          0.3,
          0.1,
          0.2,
          0.3,
          0.1,
          0.2,
          0.3,
          0.1,
          0.2,
          0.3,
          0.1,
          0.2,
          0.3)
  ref <- seq(1, 9, 1)
  min_speciation_gen <- c(1, 1, 1, 5, 5, 5, 10, 10, 10,
                          1, 1, 1, 5, 5, 5, 10, 10, 10)
  max_speciation_gen <- c(10, 10, 10, 10, 10, 10, 10, 10, 10,
                          20, 20, 20, 20, 20, 20, 20, 20, 20)
  expect_equal(sr, crefs$speciation_rate)
  expect_equal(min_speciation_gen, crefs$min_speciation_gen)
  expect_equal(max_speciation_gen, crefs$max_speciation_gen)
  expect_equal(30, tmp$getSpeciesRichness(1))
  expect_equal(42, tmp$getSpeciesRichness(2))
  expect_equal(58, tmp$getSpeciesRichness(3))
  expect_equal(21, tmp$getSpeciesRichness(4))
  expect_equal(24, tmp$getSpeciesRichness(5))
  expect_equal(24, tmp$getSpeciesRichness(6))
  
})

context("Basic spatial protracted coalescence simulations")
test_that("Basic simulation produces 100 species with low speciation generation",
          {
            tmp <- ProtractedSpatialTreeSimulation$new()
            tmp$setSimulationParameters(
              task = 203,
              seed = 3,
              min_speciation_rate = 0.000001,
              deme = 1,
              uses_logging = FALSE,
              min_speciation_gen = 0.0,
              max_speciation_gen = 0.01,
              sigma = 1,
              fine_map_file = "null",
              fine_map_x_size = 10,
              fine_map_y_size = 10
            )
            expect_equal(TRUE, tmp$runSimulation())
            tmp$applySpeciationRates(
              speciation_rates = 0.000001,
              min_speciation_gens = 0.0,
              max_speciation_gens = 0.01
            )
            expect_equal(100, tmp$getSpeciesRichness())
          })



test_that("Basic simulation produces 1 species with high speciation generation",
          {
            tmp <- ProtractedSpatialTreeSimulation$new()
            tmp$setSimulationParameters(
              task = 202,
              seed = 3,
              min_speciation_rate = 0.9999,
              deme = 1,
              sigma = 1,
              uses_logging = FALSE,
              min_speciation_gen = 1000000.0,
              max_speciation_gen = 11000000.0,
              fine_map_file = "null",
              fine_map_x_size = 10,
              fine_map_y_size = 10
            )
            expect_equal(TRUE, tmp$runSimulation())
            tmp$applySpeciationRates(
              speciation_rates = 0.9999,
              min_speciation_gens = c(1000000.0),
              max_speciation_gens = c(11000000.0)
            )
            expect_equal(1, tmp$getSpeciesRichness())
          })

test_that("Basic simulation produces 100 species with low speciation generation",
          {
            tmp <- ProtractedSpatialTreeSimulation$new()
            tmp$setSimulationParameters(
              task = 203,
              seed = 3,
              min_speciation_rate = 0.000001,
              deme = 1,
              sigma = 1,
              uses_logging = FALSE,
              min_speciation_gen = 0.0,
              max_speciation_gen = 0.01,
              fine_map_file = "null",
              fine_map_x_size = 10,
              fine_map_y_size = 10
            )
            expect_equal(TRUE, tmp$runSimulation())
            tmp$applySpeciationRates(
              speciation_rates = 0.000001,
              min_speciation_gens = 0.0,
              max_speciation_gens = 0.01
            )
            expect_equal(100, tmp$getSpeciesRichness())
          })
context("Basic spatial coalescence simulations with a metacommunity")
test_that("Metacommunity application works as intended with size of 1." , {
  tmp <- SpatialTreeSimulation$new()
  tmp$setSimulationParameters(
    task = 301,
    seed = 3,
    min_speciation_rate = 0.000001,
    deme = 1,
    sigma = 1,
    uses_logging = FALSE,
    fine_map_file = "null",
    fine_map_x_size = 10,
    fine_map_y_size = 10
  )
  expect_equal(TRUE, tmp$runSimulation())
  tmp$applySpeciationRates(
    speciation_rates = 0.000001,
    metacommunity_option = "simulated",
    metacommunity_size = 1,
    metacommunity_speciation_rate = 0.99999999,
    metacommunity_external_reference = 0
  )
  expect_equal(1, tmp$getSpeciesRichness())
})

test_that("Metacommunity application works as intended with large size." , {
  tmp <- SpatialTreeSimulation$new()
  tmp$setSimulationParameters(
    task = 302,
    seed = 3,
    min_speciation_rate = 0.000001,
    deme = 1,
    sigma = 1,
    uses_logging = FALSE,
    fine_map_file = "null",
    fine_map_x_size = 10,
    fine_map_y_size = 10
  )
  expect_equal(TRUE, tmp$runSimulation())
  tmp$applySpeciationRates(
    speciation_rates = 0.99,
    metacommunity_option = "simulated",
    metacommunity_size = 1000000,
    metacommunity_speciation_rate = 0.00001,
    metacommunity_external_reference = 0
  )
  expect_equal(23, tmp$getSpeciesRichness())
})

test_that("Metacommunity application works as intended with multiple application." ,
          {
            tmp <- SpatialTreeSimulation$new()
            tmp$setSimulationParameters(
              task = 303,
              seed = 3,
              min_speciation_rate = 0.999999,
              deme = 1,
              sigma = 1,
              uses_logging = FALSE,
              fine_map_file = "null",
              fine_map_x_size = 10,
              fine_map_y_size = 10
            )
            expect_equal(TRUE, tmp$runSimulation())
            tmp$applySpeciationRates(
              speciation_rates = 0.999999,
              metacommunity_option = c("simulated", "simulated", "simulated", "analytical"),
              metacommunity_size = c(1, 1, 1000000, 1000000),
              metacommunity_speciation_rate = c(0.999999, 0.0000000001, 0.999999, 0.0000000001),
              metacommunity_external_reference = c(0, 0, 0, 0)
            )
            tmp$output()
            expected_metacommunity = data.frame(matrix(
              c(
                1,
                0.999999,
                1,
                "simulated",
                0,
                2,
                0.0000000001,
                1,
                "simulated",
                0,
                3,
                0.999999,
                1000000,
                "simulated",
                0,
                4,
                0.0000000001,
                1000000,
                "analytical",
                0
              ),
              nrow = 4,
              ncol = 5,
              byrow = TRUE
            ), stringsAsFactors = FALSE)
            names(expected_metacommunity) <-
              c(
                "reference",
                "speciation_rate",
                "metacommunity_size",
                "option",
                "external_reference"
              )
            expected_metacommunity$reference <-
              as.numeric(expected_metacommunity$reference)
            expected_metacommunity$speciation_rate <- as.numeric(expected_metacommunity$speciation_rate)
            expected_metacommunity$metacommunity_size <-
              as.numeric(expected_metacommunity$metacommunity_size)
            expected_metacommunity$option <-
              as.character(expected_metacommunity$option)
            expected_metacommunity$external_reference <-
              as.numeric(expected_metacommunity$external_reference)
            expect_equal(TRUE,
                         all.equal(expected_metacommunity, tmp$getMetacommunityReferences()))
            expect_equal(1, tmp$getSpeciesRichness(community_reference = 1))
            expect_equal(1, tmp$getSpeciesRichness(community_reference = 2))
            expect_equal(100, tmp$getSpeciesRichness(community_reference = 3))
            expect_equal(7, tmp$getSpeciesRichness(community_reference = 4))
          })

unlink("output", recursive = TRUE)
unlink("default", recursive = TRUE)
