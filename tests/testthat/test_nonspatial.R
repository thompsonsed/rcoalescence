# Tests the non-spatial simulations

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
  expect_equal(26, tmp$getSpeciesRichness())
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
unlink(file.path("output"), recursive = TRUE)
unlink(file.path("default"), recursive = TRUE)