.onLoad <- function(libname, pkgname) {
  modules <- paste0("stan_fit4", names(stanmodels), "_mod")
  for (m in modules) loadModule(m, what = TRUE)
}

.onAttach <- function(...) {
    pkg_desc <- utils::packageDescription("paramedic")
    packageStartupMessage(paste0(
        "paramedic version ", pkg_desc$Version,
        ": ", pkg_desc$Title
    ))
    packageStartupMessage(paste0(
        "Package created on ", utils::packageDate("paramedic")
    ))
}