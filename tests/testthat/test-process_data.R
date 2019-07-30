## test that processing data works!
## set up the data
library("testthat")
library("paramedic")
data(simple_example_data)
q <- 3
q_obs <- 2

## process the data
test_that("processing works", {
  processed_data <- paramedic::process_data(full_data = simple_example_data, br_inds = 1:q,
                                 qpcr_inds = (q + 1):(q + q_obs),
                                 pcr_plus_br_inds = 1:q_obs,
                                 regex_thr = "", regex_cps = "_cps", llod = 0,
                                 m_min = 1000, div_num = 1)
  expect_equal(dim(processed_data$br)[2], q)
  expect_equal(dim(processed_data$qpcr)[2], q_obs)
})
