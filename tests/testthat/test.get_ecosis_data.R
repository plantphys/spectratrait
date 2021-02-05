context("Grabbing data via the EcoSIS API works")

test_that("Downloading data from EcoSIS doesnt throw an error", {
  ecosis_id <- "960dbb0c-144e-4563-8117-9e23d14f4aa9"
  dat_raw <- get_ecosis_data(ecosis_id = ecosis_id)
  expect_true(!is.null(dat_raw))
})