# Tests that coalescence simulations run as expected

context("Basic spatial coalescence simulations")
test_that("Basic simulation on a null landscape completes", {
  tmp <- new("SpatialTree")
  tmp$setKeyParameters(0, 0, "default", 10, 1, c(0.0))
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


context("Basic non-spatial coalescence simulations")
test_that("Basic simulation with deme of 100 completes", {
  tmp <- Tree$new()
  tmp$setKeyParameters(101, 0, "default", 1, 1, c(0.0), deme=100)
  tmp$setSpeciationParameters(0.1)
  tmp$setup()
  expect_equal(TRUE, tmp$runSimulation())
  tmp$applySpeciationRates(speciation_rates=0.1)
  expect_equal(24, tmp$getSpeciesRichness())
})

test_that("Basic simulation produces 100 species with maximum speciation rate", {
  tmp <- Tree$new()
  tmp$setSimulationParameters(task=102, seed=1, speciation_rate=0.99999, deme=100,
                              uses_logging = FALSE)
  expect_equal(TRUE, tmp$runSimulation())
  tmp$applySpeciationRates(speciation_rates=0.999999)
  expect_equal(100, tmp$getSpeciesRichness())
})

test_that("Basic simulation produces 1 species with minimum speciation rate", {
  tmp <- Tree$new()
  tmp$setSimulationParameters(task=103, seed=2, speciation_rate=0.00000000001, deme=100,
                              uses_logging = FALSE)
  expect_equal(TRUE, tmp$runSimulation())
  tmp$applySpeciationRates(speciation_rates=0.00000000001)
  expect_equal(1, tmp$getSpeciesRichness())
})



context("Basic non-spatial protracted coalescence simulations")
test_that("Basic simulation with deme of 100 completes", {
  tmp <- ProtractedTree$new()
  tmp$setKeyParameters(201, 1, "output", 1, 1, c(0.0), deme=100)
  tmp$setSpeciationParameters(0.01, min_speciation_gen = 2.0, max_speciation_gen=4000.0)
  tmp$setup()
  expect_equal(TRUE, tmp$runSimulation())
  tmp$applySpeciationRates(speciation_rates=c(0.01, 0.02),
                           min_speciation_gens=c(2.0), max_speciation_gens=c(4000.0))
  tmp$output()
  expect_equal(6, tmp$getSpeciesRichness())
  expect_equal(4, tmp$getSpeciesRichness(1))
  expect_equal(6, tmp$getSpeciesRichness(2))
})

test_that("Basic simulation produces 1 species with high speciation generation", {
  tmp <- ProtractedTree$new()
  tmp$setSimulationParameters(task=202, seed=3, speciation_rate=0.9999, deme=100,
                              uses_logging = FALSE, min_speciation_gens = c(1000000.0),
                              max_speciation_gens = c(11000000.0))
  expect_equal(TRUE, tmp$runSimulation())
  tmp$applySpeciationRates(speciation_rates=0.9999, min_speciation_gens = c(1000000.0),
                           max_speciation_gens = c(11000000.0))
  expect_equal(1, tmp$getSpeciesRichness())
})

test_that("Basic simulation produces 100 species with low speciation generation", {
  tmp <- ProtractedTree$new()
  tmp$setSimulationParameters(task=203, seed=3, speciation_rate=0.000001, deme=100,
                              uses_logging = FALSE, min_speciation_gens = 0.0,
                              max_speciation_gens = 0.01)
  expect_equal(TRUE, tmp$runSimulation())
  tmp$applySpeciationRates(speciation_rates=0.000001, min_speciation_gens = 0.0,
                           max_speciation_gens = 0.01)
  expect_equal(100, tmp$getSpeciesRichness())
})

test_that("Basic simulation produces 100 species with maximum speciation rate", {
  tmp <- ProtractedTree$new()
  tmp$setSimulationParameters(task=204, seed=1, speciation_rate=0.1, deme=100,
                              uses_logging = FALSE, min_speciation_gens=1.0, max_speciation_gens=30.0)
  expect_equal(TRUE, tmp$runSimulation())
  tmp2 <- ProtractedTree$new()
  tmp2$setSimulationParameters(task=205, seed=1, speciation_rate=0.1, deme=100,
                              uses_logging = FALSE, min_speciation_gens=1.0, max_speciation_gens=10.0)
  expect_equal(TRUE, tmp2$runSimulation())
  tmp3 <- ProtractedTree$new()
  tmp3$setSimulationParameters(task=206, seed=1, speciation_rate=0.1, deme=100,
                               uses_logging = FALSE, min_speciation_gens=5.0, max_speciation_gens=10.0)
  expect_equal(TRUE, tmp3$runSimulation())
  tmp$applySpeciationRates(speciation_rates=0.1, min_speciation_gens=1.0, max_speciation_gens=30.0)
  r1 <- tmp$getSpeciesRichness()
  tmp2$applySpeciationRates(speciation_rates=0.1, min_speciation_gens=1.0, max_speciation_gens=10.0)
  r2 <- tmp2$getSpeciesRichness()
  tmp3$applySpeciationRates(speciation_rates=0.1, min_speciation_gens=5.0, max_speciation_gens=10.0)
  r3 <- tmp3$getSpeciesRichness()
  expect_equal(TRUE, r2 > r1)
  expect_equal(TRUE, r2 > r3)
  expect_equal(28, r2)
  expect_equal(25, r1)
  expect_equal(20, r3)
})


test_that("Can apply multiple protracted speciation parameters", {
  tmp <- ProtractedTree$new()
  tmp$setSimulationParameters(task=204, seed=1, speciation_rate=0.1, deme=100,
                              uses_logging = FALSE, min_speciation_gens=10.0, max_speciation_gens=30.0)
  expect_equal(TRUE, tmp$runSimulation())
  expect_error(tmp$applySpeciationRates(speciation_rates=c(0.1, 0.2, 0.3),
                                        min_speciation_gens=c(1.0, 5.0, 10.0), 
                                        max_speciation_gens=c(10.0, 10.0)))
  tmp$applySpeciationRates(speciation_rates=c(0.1, 0.2, 0.3),
                           min_speciation_gens=c(1.0, 5.0, 10.0), 
                           max_speciation_gens=c(10.0, 10.0, 20.0))
  tmp$output()
  crefs <- tmp$getCommunityReferences()
  sr <- c(0.1, 0.2, 0.3, 0.1, 0.2, 0.3, 0.1, 0.2, 0.3)
  ref <- seq(1, 9, 1)
  min_speciation_gen <- c(1, 1, 1, 5, 5, 5, 10, 10, 10)
  max_speciation_gen <- c(10, 10, 10, 10, 10, 10, 20, 20, 20)
  expect_equal(sr, crefs$speciation_rate)
  expect_equal(min_speciation_gen, crefs$min_speciation_gen)
  expect_equal(max_speciation_gen, crefs$max_speciation_gen)
  expect_equal(30, tmp$getSpeciesRichness(1))
  expect_equal(42, tmp$getSpeciesRichness(2))
  expect_equal(58, tmp$getSpeciesRichness(3))
  expect_equal(21, tmp$getSpeciesRichness(4))
  expect_equal(24, tmp$getSpeciesRichness(5))
  expect_equal(24, tmp$getSpeciesRichness(6))
  expect_equal(13, tmp$getSpeciesRichness(7))
  expect_equal(15, tmp$getSpeciesRichness(8))
  expect_equal(15, tmp$getSpeciesRichness(9))
  
})

context("Basic spatial protracted coalescence simulations")
test_that("Basic simulation produces 100 species with low speciation generation", {
  tmp <- ProtractedSpatialTree$new()
  tmp$setSimulationParameters(task=203, seed=3, speciation_rate=0.000001, deme=1,
                              uses_logging = FALSE, min_speciation_gens = 0.0,
                              max_speciation_gens = 0.01, sigma=1, fine_map_file="null",
                              fine_map_x_size = 10, fine_map_y_size = 10)
  expect_equal(TRUE, tmp$runSimulation())
  tmp$applySpeciationRates(speciation_rates=0.000001, min_speciation_gens = 0.0,
                           max_speciation_gens = 0.01)
  expect_equal(100, tmp$getSpeciesRichness())
})



test_that("Basic simulation produces 1 species with high speciation generation", {
  tmp <- ProtractedSpatialTree$new()
  tmp$setSimulationParameters(task=202, seed=3, speciation_rate=0.9999, deme=1, sigma=1,
                              uses_logging = FALSE, min_speciation_gens = c(1000000.0),
                              max_speciation_gens = c(11000000.0), fine_map_file="null",
                              fine_map_x_size = 10, fine_map_y_size = 10)
  expect_equal(TRUE, tmp$runSimulation())
  tmp$applySpeciationRates(speciation_rates=0.9999, min_speciation_gens = c(1000000.0),
                           max_speciation_gens = c(11000000.0))
  expect_equal(1, tmp$getSpeciesRichness())
})

test_that("Basic simulation produces 100 species with low speciation generation", {
  tmp <- ProtractedSpatialTree$new()
  tmp$setSimulationParameters(task=203, seed=3, speciation_rate=0.000001, deme=1, sigma=1,
                              uses_logging = FALSE, min_speciation_gens = 0.0,
                              max_speciation_gens = 0.01, fine_map_file="null",
                              fine_map_x_size = 10, fine_map_y_size = 10)
  expect_equal(TRUE, tmp$runSimulation())
  tmp$applySpeciationRates(speciation_rates=0.000001, min_speciation_gens = 0.0,
                           max_speciation_gens = 0.01)
  expect_equal(100, tmp$getSpeciesRichness())
})


unlink("output", recursive=TRUE)

