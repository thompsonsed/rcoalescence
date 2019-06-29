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
    fine_map_file = "sample/example_fine.tif",
    fine_map_x_size = 13,
    fine_map_y_size = 13,
    coarse_map_file = "sample/example_coarse.tif",
    coarse_map_x_size = 35,
    coarse_map_y_size = 41,
    coarse_map_x_offset = 11,
    coarse_map_y_offset = 14,
    sample_mask_file = "sample/example_mask.tif",
    sample_x_size = 13,
    sample_y_size = 13,
    uses_logging = FALSE
  )
  expect_equal(TRUE, tmp$runSimulation())
  tmp$applySpeciationRates(speciation_rates = 0.5)
  expect_equal(1160, tmp$getSpeciesRichness())
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
    fine_map_file = "sample/example_fine.tif",
    coarse_map_file = "sample/example_coarse.tif",
    sample_mask_file = "sample/example_mask.tif",
    uses_logging = FALSE
  )
  expect_equal(TRUE, tmp$runSimulation())
  tmp$applySpeciationRates(speciation_rates = c(0.5))
  expect_equal(1160, tmp$getSpeciesRichness())
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
        coarse_map_file = "sample/example_coarse.tif",
        coarse_map_x = 35L,
        coarse_map_y = 41L,
        coarse_map_x_offset = 11L,
        coarse_map_y_offset = 14L,
        coarse_map_scale = 1,
        fine_map_file = "sample/example_fine.tif",
        fine_map_x = 13L,
        fine_map_y = 13L,
        fine_map_x_offset = 0L,
        fine_map_y_offset = 0L,
        sample_file = "sample/example_mask.tif",
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
      row.names = c(NA,-1L)
    )
  expect_equal(TRUE, all.equal(sim_params, expected_params))
  
})
unlink(file.path("output"), recursive = TRUE)
unlink(file.path("default"), recursive = TRUE)