# Tests that coalescence simulations run as expected
context("Checking biodiversity metrics simulations")
test_that("Biodiversity metrics correctly stored in output database", {
  output_file <- file.path("output", "data_1_10.db")
  if (file.exists(output_file)) {
    file.remove(output_file)
  }
  tmp <- SpatialTreeSimulation$new()
  tmp$setSimulationParameters(
    seed = 10,
    task = 1,
    output_directory = "output",
    min_speciation_rate = 0.5,
    sigma = 2 * (2**0.5),
    deme = 1,
    deme_sample = 0.1,
    fine_map_file = "sample/example_fine.tif",
    coarse_map_file = "sample/example_coarse.tif",
    sample_mask_file = "sample/example_mask.tif",
    uses_logging = FALSE
  )
  expect_equal(TRUE, tmp$runSimulation())
  tmp$applySpeciationRates(
    speciation_rates = c(0.5, 0.7),
    use_spatial = TRUE,
    record_ages = TRUE
  )
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
    c(
      "reference",
      "speciation_rate",
      "time",
      "fragments",
      "metacommunity_reference"
    )
  expect_equal(
    TRUE,
    all.equal(tmp$getCommunityReferences(), community_references)
  )
  # Check species abundances
  expected_abundances <-
    data.frame(matrix(c(0, 0, 1, 1, 2, 1, 3, 1, 4, 1), ncol = 2, byrow = TRUE))
  names(expected_abundances) <- c("species_id", "no_individuals")
  actual_abundances <- head(tmp$getSpeciesAbundances(1), 5)
  expect_equal(TRUE, all.equal(expected_abundances, actual_abundances))
  expected_abundances <-
    data.frame(matrix(
      c(1163, 2, 1164, 1, 1165, 1, 1166, 1, 1167, 1),
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
      c(1, 0, 10, 2, 0, 10, 3, 0, 10, 4, 0, 10, 1147, 0, 10),
      ncol = 3,
      byrow = TRUE
    ))
  names(expected_locations) <- c("species_id", "x", "y")
  actual_locations <- head(tmp$getSpeciesLocations(1), 5)
  expect_equal(TRUE, all.equal(expected_locations, actual_locations))
  expected_locations <-
    data.frame(matrix(
      c(1147, 12, 12, 1148, 12, 12, 1149, 12, 12, 1150, 12, 12, 1151, 12, 12),
      ncol = 3,
      byrow = TRUE
    ))
  names(expected_locations) <- c("species_id", "x", "y")
  actual_locations <-
    data.frame(tail(tmp$getSpeciesLocations(2), 5), row.names = NULL)
  expect_equal(TRUE, isTRUE(all.equal(expected_locations, actual_locations)))
  # Check species richness
  expect_equal(1157, tmp$getSpeciesRichness(1))
  expect_equal(1167, tmp$getSpeciesRichness(2))
  expected_ages <- read.csv("sample/species_ages_example_results.csv") %>%
    select(-X)
  expect_equal(expected_ages, tmp$getAllSpeciesAges())
  expected_df <- read.csv("sample/species_ages_example_results_single.csv") %>% select(-X)
  expect_equal(expected_df, tmp$getSpeciesAges())
})

context("More complex spatial coalescence simulations")
test_that("Simulation with a single historical maps works as intended.", {
  tmp <- SpatialTreeSimulation$new()
  tmp$setSimulationParameters(
    seed = 10,
    task = 1,
    output_directory = "output",
    min_speciation_rate = 0.1,
    sigma = 2 * (2**0.5),
    deme = 1,
    deme_sample = 0.1,
    fine_map_file = "sample/example_fine.tif",
    coarse_map_file = "sample/example_coarse.tif",
    sample_mask_file = "sample/example_mask.tif",
    uses_logging = FALSE,
    landscape_type = "closed",
    partial_setup = TRUE
  )
  tmp$addHistoricalMap(
    historical_fine_map = "sample/example_historical_fine.tif",
    historical_coarse_map = "sample/example_coarse.tif",
    gen_since_historical = 1,
    habitat_change_rate = 0.5
  )
  # tmp$addHistoricalMap(historical_fine_map="sample/example_fine.tif",
  #                      historical_coarse_map="sample/example_coarse.tif",
  #                      gen_since_historical=10, habitat_change_rate=0.5)
  # tmp$addHistoricalMap(historical_fine_map="sample/example_fine.tif",
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
    sigma = 2 * (2**0.5),
    deme = 1,
    deme_sample = 0.1,
    fine_map_file = "sample/example_fine.tif",
    coarse_map_file = "sample/example_coarse.tif",
    sample_mask_file = "sample/example_mask.tif",
    uses_logging = FALSE,
    partial_setup = TRUE
  )
  tmp$addHistoricalMap(
    historical_fine_map = "sample/example_historical_fine.tif",
    historical_coarse_map = "sample/example_coarse.tif",
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
    historical_fine_map = "sample/example_fine.tif",
    historical_coarse_map = "sample/example_coarse.tif",
    gen_since_historical = 20,
    habitat_change_rate = 1.0
  )
  expect_equal(TRUE, tmp$runSimulation())
  tmp$applySpeciationRates(speciation_rates = c(0.5))
  expect_equal(1157, tmp$getSpeciesRichness())
})
#
test_that("Simulation with multiple sampling times.", {
  tmp <- SpatialTreeSimulation$new()
  tmp$setSimulationParameters(
    seed = 11,
    task = 1,
    output_directory = "output",
    min_speciation_rate = 0.5,
    sigma = 2 * (2**0.5),
    deme = 1,
    deme_sample = 0.1,
    fine_map_file = "sample/example_fine.tif",
    coarse_map_file = "sample/example_coarse.tif",
    sample_mask_file = "sample/example_mask.tif",
    uses_logging = FALSE,
    times_list = c(0.0, 1.0, 10.0, 20.0)
  )
  expect_equal(TRUE, tmp$runSimulation())
  tmp$applySpeciationRates(
    speciation_rates = c(0.5),
    times_list = c(0.0, 1.0, 10.0, 20.0)
  )
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
    c(
      "reference",
      "speciation_rate",
      "time",
      "fragments",
      "metacommunity_reference"
    )
  expect_equal(
    TRUE,
    all.equal(community_references, tmp$getCommunityReferences())
  )
  expect_equal(1153, tmp$getSpeciesRichness(1))
  expect_equal(1164, tmp$getSpeciesRichness(2))
  expect_equal(1167, tmp$getSpeciesRichness(3))
  expect_equal(1167, tmp$getSpeciesRichness(4))
})

unlink(file.path("output"), recursive = TRUE)
unlink(file.path("default"), recursive = TRUE)
