# paramedic 0.1.3

## Major changes

* added a second option for obtaining predictions on new data using a previously fitted model: the options are now `predict_paramedic` and `posterior_predict`. The former uses posterior draws and a Stan algorithm to generate predictions; the latter uses posterior draws and internal R functions to generate predictions. Both are expected to perform similarly.

## Minor changes

* Updated return type of `run_paramedic` to class `paramedic`; added associated `print`, `summary` functions

# paramedic 0.1.2

## Major changes

* added ability to obtain predictions on new data using a previously fitted model using `predict_paramedic`

## Minor changes

None

# paramedic 0.1.1

## Major changes

* added ability to run analysis adjusted for batch effects by setting argument `k` to be nonzero in `run_paramedic` (the default is zero for data collected in a single experiment)
* added new stan code for batch analysis

## Minor changes

* reorganized stan code into directories, since much of the code was shared between models; the new organization uses `#include` statements, and the shared code is stored in `data/` (sets up data objects), `parameters/` (sets up model parameters), `tparameters/` (sets up transformed parameters), and `model/` (sets up hierarchical model for shared parameters)

# paramedic 0.1.0

## Major changes

* removed function `run_paramedic_centered`; the desired behavior can now be obtained by using `run_paramedic` with the argument `"centered = TRUE"`
* removed function `no_efficiency`; the desired behavior can now be obtained by using `run_paramedic` with arguments `"alpha_sigma = 0"` and `"kappa_sigma = 0"`
* removed several stan models by collapsing functionality into `variable_efficiency.stan` and `variable_efficiency_centered.stan`

## Minor changes

* updated all tests to reflect major changes
* added vignette discussing hierarchical model specification
* added a `pkgdown` site

# paramedic 0.0.3.900

## Major changes

None

## Minor changes

* new internal functions to reduce replicating code within key functions
* new `stan` model for specifying a negative binomial distribution on V

# paramedic 0.0.3

## Major changes

None

## Minor changes

* update internal `Stan` models to use log densities and unconstrained parameters, then back-transform after sampling

# paramedic 0.0.2

## Major changes

* allow an option to adjust for covariates

## Minor changes

None

# paramedic 0.0.1

## Major changes

None

## Minor changes

* allow `run_paramedic` with `q` = `q_obs`

# paramedic 0.0.0.9000

This is the initial package release!
