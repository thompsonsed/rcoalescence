context("Basic non-spatial coalescence simulations with a metacommunity")

test_that("Metacommunity application works as intended with non-spatial local communities." , {
  tmp <- TreeSimulation$new()
  tmp$setSimulationParameters(
    task = 300,
    seed = 10,
    min_speciation_rate = 0.00001,
    deme = 100,
    uses_logging = FALSE
  )
  expect_equal(TRUE, tmp$runSimulation())
  tmp$setSpeciationParameters(speciation_rates = c(0.00001, 0.99),
                              metacommunity_option = c("simulated", "analytical"),
                              metacommunity_size = c(100000, 1000000),
                              metacommunity_speciation_rate = c(0.00001, 0.99))
  tmp$applySpeciationRates()
  expect_equal(100, tmp$getSpeciesRichness())
  tmp$output()
  expect_equal(1, tmp$getSpeciesRichness(1))
  expect_equal(2, tmp$getSpeciesRichness(2))
  expect_equal(1, tmp$getSpeciesRichness(3))
  expect_equal(100, tmp$getSpeciesRichness(4))
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
  expect_equal(32, tmp$getSpeciesRichness())
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
              fine_map_y_size = 10,
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
            expected_metacommunity$speciation_rate <-
              as.numeric(expected_metacommunity$speciation_rate)
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
            expect_equal(1, tmp$getSpeciesRichness(community_reference = 4))
          })

unlink(file.path("output"), recursive = TRUE)
unlink(file.path("default"), recursive = TRUE)