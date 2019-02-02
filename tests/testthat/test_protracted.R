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

context("Protracted simulation sanity checks")
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

context("More complex spatial protracted simulations")
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

test_that("More complex spatial example using protracted speciation",
          {
            tmp <- ProtractedSpatialTreeSimulation$new()
            tmp$setSimulationParameters(
              seed = 1,
              task = 204,
              output_directory = "output",
              min_speciation_rate = 0.01,
              sigma = 4,
              deme = 10,
              deme_sample = 0.001,
              fine_map_file = "sample/example_fine.tif",
              coarse_map_file = "sample/example_coarse.tif",
              sample_mask_file = "sample/example_mask.tif",
              landscape_type = "tiled_fine",
              uses_logging = FALSE,
              min_speciation_gen = 10,
              max_speciation_gen = 100000,
              times_list = c(0.0, 1.0, 10.0, 20.0)
            )
            tmp$addHistoricalMap(
              historical_fine_map = "sample/example_historical_fine.tif",
              historical_coarse_map = "sample/example_coarse.tif",
              gen_since_historical = 5,
              habitat_change_rate = 0.0
            )
            tmp$addHistoricalMap(
              historical_fine_map = "sample/example_fine.tif",
              historical_coarse_map = "sample/example_coarse.tif",
              gen_since_historical = 10,
              habitat_change_rate = 0.0
            )
            tmp$addHistoricalMap(
              historical_fine_map = "sample/example_historical_fine.tif",
              historical_coarse_map = "sample/example_coarse.tif",
              gen_since_historical = 20,
              habitat_change_rate = 0.0
            )
            expect_equal(TRUE, tmp$runSimulation())
            tmp$applySpeciationRates(
              speciation_rates = c(0.01, 0.02),
              times_list = c(0.0, 1.0, 10.0, 20.0),
              min_speciation_gens = c(0, 0, 10, 10),
              max_speciation_gens = c(1000, 10000, 10000, 100000)
            )
            tmp$output()
            community_references <-
              structure(
                list(
                  reference = 1:32,
                  speciation_rate = c(
                    0.01,
                    0.01,
                    0.01,
                    0.01,
                    0.02,
                    0.02,
                    0.02,
                    0.02,
                    0.01,
                    0.01,
                    0.01,
                    0.01,
                    0.02,
                    0.02,
                    0.02,
                    0.02,
                    0.01,
                    0.01,
                    0.01,
                    0.01,
                    0.02,
                    0.02,
                    0.02,
                    0.02,
                    0.01,
                    0.01,
                    0.01,
                    0.01,
                    0.02,
                    0.02,
                    0.02,
                    0.02
                  ),
                  time = c(
                    0,
                    1,
                    10,
                    20,
                    0,
                    1,
                    10,
                    20,
                    0,
                    1,
                    10,
                    20,
                    0,
                    1,
                    10,
                    20,
                    0,
                    1,
                    10,
                    20,
                    0,
                    1,
                    10,
                    20,
                    0,
                    1,
                    10,
                    20,
                    0,
                    1,
                    10,
                    20
                  ),
                  fragments = c(
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L
                  ),
                  metacommunity_reference = c(
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L,
                    0L
                  ),
                  min_speciation_gen = c(
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    10,
                    10,
                    10,
                    10,
                    10,
                    10,
                    10,
                    10,
                    10,
                    10,
                    10,
                    10,
                    10,
                    10,
                    10,
                    10
                  ),
                  max_speciation_gen = c(
                    1000,
                    1000,
                    1000,
                    1000,
                    1000,
                    1000,
                    1000,
                    1000,
                    10000,
                    10000,
                    10000,
                    10000,
                    10000,
                    10000,
                    10000,
                    10000,
                    10000,
                    10000,
                    10000,
                    10000,
                    10000,
                    10000,
                    10000,
                    10000,
                    1e+05,
                    1e+05,
                    1e+05,
                    1e+05,
                    1e+05,
                    1e+05,
                    1e+05,
                    1e+05
                  )
                ),
                .Names = c(
                  "reference",
                  "speciation_rate",
                  "time",
                  "fragments",
                  "metacommunity_reference",
                  "min_speciation_gen",
                  "max_speciation_gen"), class = "data.frame", row.names = c(NA, 
                                                                                                                                                                                                                                                                                              -32L))
            expect_equal(TRUE,
                         all.equal(community_references, tmp$getCommunityReferences()))
            expect_equal(76, tmp$getSpeciesRichness(1))
            expect_equal(76, tmp$getSpeciesRichness(2))
            expect_equal(6, tmp$getSpeciesRichness(3))
            expect_equal(6, tmp$getSpeciesRichness(4))
            expect_equal(76, tmp$getSpeciesRichness(5))
            expect_equal(76, tmp$getSpeciesRichness(6))
            expect_equal(6, tmp$getSpeciesRichness(7))
            expect_equal(6, tmp$getSpeciesRichness(8))
            expect_equal(76, tmp$getSpeciesRichness(9))
            expect_equal(76, tmp$getSpeciesRichness(10))
            expect_equal(6, tmp$getSpeciesRichness(11))
            expect_equal(6, tmp$getSpeciesRichness(12))
            expect_equal(76, tmp$getSpeciesRichness(13))
            expect_equal(76, tmp$getSpeciesRichness(14))
            expect_equal(6, tmp$getSpeciesRichness(15))
            expect_equal(6, tmp$getSpeciesRichness(16))
            expect_equal(76, tmp$getSpeciesRichness(17))
            expect_equal(76, tmp$getSpeciesRichness(18))
            expect_equal(6, tmp$getSpeciesRichness(19))
            expect_equal(6, tmp$getSpeciesRichness(20))
            expect_equal(76, tmp$getSpeciesRichness(21))
            expect_equal(76, tmp$getSpeciesRichness(22))
            expect_equal(6, tmp$getSpeciesRichness(23))
            expect_equal(6, tmp$getSpeciesRichness(24))
            expect_equal(76, tmp$getSpeciesRichness(25))
            expect_equal(76, tmp$getSpeciesRichness(26))
            expect_equal(6, tmp$getSpeciesRichness(27))
            expect_equal(6, tmp$getSpeciesRichness(28))
            expect_equal(76, tmp$getSpeciesRichness(29))
            expect_equal(76, tmp$getSpeciesRichness(30))
            expect_equal(6, tmp$getSpeciesRichness(31))
            expect_equal(6, tmp$getSpeciesRichness(32))
            
          })

unlink(file.path("output"), recursive = TRUE)
unlink(file.path("default"), recursive = TRUE)