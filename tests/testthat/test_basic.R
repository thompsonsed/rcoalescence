# Tests that coalescence simulations run as expected

context("Basic coalescence simulations")
test_that("Basic simulation on a null landscape completes", {
  tmp <- new("SpatialTree")
  tmp$setKeyParameters(0, 0, "default", 1, 1, c(0.0))
  tmp$setSpeciationParameters(0.1)
  tmp$setPristineMapParameters()
  tmp$setMapParameters()
  tmp$setup()
  expect_equal(TRUE, tmp$runSimulation())
})

test_that("Basic simulation produces the expected number of individuals in edge cases", {
  tmp <- new("SpatialTree")
  params <- SimulationParameters
  params$task <- 12
  params$seed <- 1
  params$speciation_rate <- 0.01
  params$deme = 1
  params$sigma = 2
  params$fine_map_x_size = 10
  params$fine_map_y_size = 10
  spec_params <- SpeciationParameters
  spec_params$speciation_rates <- c(0.01)
  tmp$setParameters(params)
  expect_equal(TRUE, tmp$runSimulation())
  tmp$applySpeciationRates(spec_params)
  expect_equal(6, tmp$getSpeciesRichness())
})