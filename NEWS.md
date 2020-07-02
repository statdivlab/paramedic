# paramedic 0.1.0

## Major changes

* removed function run_paramedic_centered; the desired behavior can now be obtained by using run_paramedic with the argument "centered = TRUE"
* removed function no_efficiency; the desired behavior can now be obtained by using run_paramedic with arguments "alpha_sigma = 0" and "kappa_sigma = 0"
* removed several stan models by collapsing functionality into variable_efficiency.stan and variable_efficiency_centered.stan

## Minor changes

* updated all tests to reflect major changes
* added vignette discussing hierarchical model specification
* added a pkgdown site

# paramedic 0.0.3.900

## Major changes

None

## Minor changes

* new internal functions to reduce replicating code within key functions
* new stan model for specifying a negative binomial distribution on V

# paramedic 0.0.3

## Major changes

None

## Minor changes

* update internal Stan models to use log densities and unconstrained parameters, then back-transform after sampling

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
