# Tests that coalescence simulations run as expected

context("Basic coalescence simulations")
test_that("Basic simulation on a null landscape completes", {
  tmp <- new("SpatialTree")
  tmp$setKeyParameters(0, 0, "default", 1, 1, c(0.0))
  tmp$setSpeciationParameters(0.1)
  tmp$setHistoricalMapParameters()
  tmp$setMapParameters(fine_map_x_size = 10, fine_map_y_size=10)
  tmp$setup()
  expect_equal(TRUE, tmp$runSimulation())
})

test_that("Basic simulation produces the expected number of individuals in edge cases", {
  tmp <- SpatialTree$new()
  tmp$setSimulationParameters(task=12, seed=1, speciation_rate=0.01, deme=1, sigma=2,
                              fine_map_x_size=10, fine_map_y_size=10, uses_logging = FALSE)
  expect_equal(TRUE, tmp$runSimulation())
  tmp$applySpeciationRates(speciation_rate=0.01)
  expect_equal(6, tmp$getSpeciesRichness())
})

test_that("Running with a basic map file works as expected", {
  tmp <- SpatialTree$new()
  tmp$setSimulationParameters(seed=9, task =1, output_directory = "output", speciation_rate = 0.5,
                              sigma=2 * (2 ** 0.5), deme=1, deme_sample = 0.1, 
                              fine_map_file = "../../inst/extdata/sample/example_fine.tif", fine_map_x_size = 13,
                              fine_map_y_size = 13, coarse_map_file="../../inst/extdata/sample/example_coarse.tif",
                              coarse_map_x_size = 35, coarse_map_y_size = 41, 
                              coarse_map_x_offset = 11, coarse_map_y_offset = 14,
                              sample_mask_file = "../../inst/extdata/sample/example_mask.tif", sample_x_size = 13,
                              sample_y_size = 13, uses_logging = FALSE)
  expect_equal(TRUE, tmp$runSimulation())
  tmp$applySpeciationRates(speciation_rates = c(0.5))
  expect_equal(1163, tmp$getSpeciesRichness())
})

test_that("Simulations correctly detect map sizes", {
  tmp <- SpatialTree$new()
  tmp$setSimulationParameters(seed=9, task =1, output_directory = "output", speciation_rate = 0.5,
                              sigma=2 * (2 ** 0.5), deme=1, deme_sample = 0.1, 
                              fine_map_file = "../../inst/extdata/sample/example_fine.tif",
                              coarse_map_file="../../inst/extdata/sample/example_coarse.tif",
                              sample_mask_file = "../../inst/extdata/sample/example_mask.tif",
                              uses_logging = FALSE)
  expect_equal(TRUE, tmp$runSimulation())
  tmp$applySpeciationRates(speciation_rates = c(0.5))
  expect_equal(1163, tmp$getSpeciesRichness())
})

test_that("Output is created successfully", {
  output_file = "output/SQL_data/data_1_10.db"
  if(file.exists(output_file))
  {
    file.remove(output_file)
  }
  tmp <- SpatialTree$new()
  tmp$setSimulationParameters(seed=10, task =1, output_directory = "output", speciation_rate = 0.5,
                              sigma=2 * (2 ** 0.5), deme=1, deme_sample = 0.1,
                              fine_map_file = "../../inst/extdata/sample/example_fine.tif",
                              coarse_map_file="../../inst/extdata/sample/example_coarse.tif",
                              sample_mask_file = "../../inst/extdata/sample/example_mask.tif",
                              uses_logging = FALSE)
  expect_equal(TRUE, tmp$runSimulation())
  tmp$applySpeciationRates(speciation_rates = c(0.5, 0.7), use_spatial = TRUE)
  tmp$output()
  # Make sure an output exists
  expect_equal(TRUE, file.exists(output_file))
  # Check community references
  community_references <- data.frame(matrix(c(1, 0.5, 0.0, 0, 0, 2, 0.7, 0.0, 0, 0), nrow = 2, ncol=5,
                                            byrow=TRUE))
  names(community_references) <- c("reference", "speciation_rate", "time", "fragments", "metacommunity_reference")
  expect_equal(TRUE, all.equal(tmp$getCommunityReferences(), community_references))
  # Check species abundances
  expected_abundances <- data.frame(matrix(c(0, 0, 1, 1, 2, 1, 3, 1, 4, 1), ncol=2, byrow=TRUE))
  names(expected_abundances) <- c("species_id", "no_individuals")
  actual_abundances <- head(tmp$getSpeciesAbundances(1), 5)
  expect_equal(TRUE, all.equal(expected_abundances, actual_abundances))
  expected_abundances <- data.frame(matrix(c(1165, 2, 1166, 1, 1167, 1, 1168, 1, 1169, 1), ncol=2, byrow=TRUE))
  names(expected_abundances) <- c("species_id", "no_individuals")
  actual_abundances <- data.frame(tail(tmp$getSpeciesAbundances(2), 5), row.names = NULL)
  expect_equal(TRUE, all.equal(expected_abundances, actual_abundances))
  # Check species locations
  expected_locations <- data.frame(matrix(c(1, 0, 10, 2, 0, 10, 3, 0, 10, 4, 0, 10, 5, 0, 10), ncol=3, byrow=TRUE))
  names(expected_locations) <- c("species_id", "x", "y")
  actual_locations <- head(tmp$getSpeciesLocations(1), 5)
  expect_equal(TRUE, all.equal(expected_locations, actual_locations))
  expected_locations <- data.frame(matrix(c(1153, 12, 12, 1154, 12, 12, 1155, 12, 12, 1156, 12, 12, 1157, 12, 12),
                                          ncol=3, byrow=TRUE))
  names(expected_locations) <- c("species_id", "x", "y")
  actual_locations <- data.frame(tail(tmp$getSpeciesLocations(2), 5), row.names=NULL)
  expect_equal(TRUE, all.equal(expected_locations, actual_locations))
  # Check species richness
  expect_equal(1161, tmp$getSpeciesRichness(1))
  expect_equal(1169, tmp$getSpeciesRichness(2))
  
  
  
})
# 
test_that("Simulation with a single historical maps works as intended.", {
  tmp <- SpatialTree$new()
  tmp$setSimulationParameters(seed=10, task =1, output_directory = "output", speciation_rate = 0.1,
                              sigma=2 * (2 ** 0.5), deme=1, deme_sample = 0.1,
                              fine_map_file = "../../inst/extdata/sample/example_fine.tif",
                              coarse_map_file="../../inst/extdata/sample/example_coarse.tif",
                              sample_mask_file = "../../inst/extdata/sample/example_mask.tif",
                              uses_logging = FALSE, landscape_type = "closed")
  tmp$addHistoricalMap(historical_fine_map="../../inst/extdata/sample/example_historical_fine.tif",
  historical_coarse_map="../../inst/extdata/sample/example_coarse.tif", gen_since_historical=1,
  habitat_change_rate=0.5)
  # tmp$addHistoricalMap(historical_fine_map="../../inst/extdata/sample/example_fine.tif",
  #                      historical_coarse_map="../../inst/extdata/sample/example_coarse.tif", 
  #                      gen_since_historical=10, habitat_change_rate=0.5)
  # tmp$addHistoricalMap(historical_fine_map="../../inst/extdata/sample/example_fine.tif",
                       # historical_coarse_map="null", gen_since_historical=4,
                       # habitat_change_rate=1.0)
  expect_equal(TRUE, tmp$runSimulation())
  tmp$applySpeciationRates(speciation_rates = c(0.1))
  expect_equal(1102, tmp$getSpeciesRichness())
})
# 
test_that("Simulation with multiple historical maps works as intended.", {
  tmp <- SpatialTree$new()
  tmp$setSimulationParameters(seed=10, task =1, output_directory = "output", speciation_rate = 0.5,
                              sigma=2 * (2 ** 0.5), deme=1, deme_sample = 0.1,
                              fine_map_file = "../../inst/extdata/sample/example_fine.tif",
                              coarse_map_file="../../inst/extdata/sample/example_coarse.tif",
                              sample_mask_file = "../../inst/extdata/sample/example_mask.tif",
                              uses_logging = FALSE)
  tmp$addHistoricalMap(historical_fine_map="../../inst/extdata/sample/example_historical_fine.tif",
                       historical_coarse_map="../../inst/extdata/sample/example_coarse.tif",
                       gen_since_historical=1, habitat_change_rate=0.5)
  tmp$addHistoricalMap(historical_fine_map="null", historical_coarse_map="null",
                       gen_since_historical=10, habitat_change_rate=1.0)
  tmp$addHistoricalMap(historical_fine_map="../../inst/extdata/sample/example_fine.tif",
                       historical_coarse_map="../../inst/extdata/sample/example_coarse.tif",
                       gen_since_historical=20, habitat_change_rate=1.0)
  expect_equal(TRUE, tmp$runSimulation())
  tmp$applySpeciationRates(speciation_rates = c(0.5))
  expect_equal(1158, tmp$getSpeciesRichness())
})
# 
test_that("Simulation with multiple sampling times.", {
  tmp <- SpatialTree$new()
  tmp$setSimulationParameters(seed=11, task =1, output_directory = "output", speciation_rate = 0.5,
                              sigma=2 * (2 ** 0.5), deme=1, deme_sample = 0.1,
                              fine_map_file = "../../inst/extdata/sample/example_fine.tif",
                              coarse_map_file="../../inst/extdata/sample/example_coarse.tif",
                              sample_mask_file = "../../inst/extdata/sample/example_mask.tif",
                              uses_logging = FALSE, times_list=c(0.0, 1.0, 10.0, 20.0))
  expect_equal(TRUE, tmp$runSimulation())
  tmp$applySpeciationRates(speciation_rates = c(0.5), times_list=c(0.0, 1.0, 10.0, 20.0))
  tmp$output()
  community_references <- data.frame(matrix(c(1, 0.5, 0.0, 0, 0, 2, 0.5, 1.0, 0, 0,
                                              3, 0.5, 10.0, 0, 0, 4, 0.5, 20.0, 0, 0),
                                            ncol=5, byrow=TRUE))
  names(community_references) <- c("reference", "speciation_rate", "time", "fragments",
                                   "metacommunity_reference")
  expect_equal(TRUE, all.equal(community_references, tmp$getCommunityReferences()))
  expect_equal(1158, tmp$getSpeciesRichness(1))
  expect_equal(1158, tmp$getSpeciesRichness(2))
  expect_equal(1164, tmp$getSpeciesRichness(3))
  expect_equal(1156, tmp$getSpeciesRichness(4))
})

unlink("output", recursive=TRUE)
