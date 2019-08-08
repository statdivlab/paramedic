#' Process a dataset, returning observed abolute and relative abundance matrices
#'
#' Return matrices of observed absolute and relative abundances, given a single dataset and a set of indices.
#'
#' @param full_data the full dataset, containing both absolute and relative abundances
#' @param rel_inds the indices marking the relative abundance measurements in \code{full_data}
#' @param abs_inds the indices marking the absolute abundance measurements in \code{full_data}
#' @param abs_plus_rel_inds the indices that have both absolute and relative abundance measurements
#' @param regex_thr the regular expression to identify absolute abundance threshold for each taxon-specific measurement (e.g., "_t" appended to each threshold)
#' @param regex_abs the regular expression to identify absolute abundance measurement for each taxon-specific measurement (e.g., "_cps" appended to each measurement of copies per swab)
#' @param llod the lower limit of detection for absolute abundance measurements
#' @param m_min the minimum number of reads to consider
#' @param div_num the number to divide absolute abundance measurements by (to create valid R integers for use in Stan)
#'
#' @return a list containing two matrices: one with the observed absolute abundances, and the other with the observed relative abundances.
#'
#' @examples
#' ## load the package, read in example data
#' library("paramedic")
#' data(simple_example_data)
#'
#' ## process the example data
#' q <- 3
#' q_obs <- 2
#' processed_data <- process_data(full_data = simple_example_data, rel_inds = 1:q,
#' abs_inds = (q + 1):(q + 1 + q_obs),
#' abs_plus_rel_inds = 1:q_obs,
#' regex_thr = "", regex_abs = "", llod = 0,
#' m_min = 1000, div_num = 1)
#'
#' @export
process_data <- function(full_data, rel_inds, abs_inds, abs_plus_rel_inds, regex_thr = "_t", regex_abs = "_cps", llod = 0, m_min = 1000, div_num = 100) {
  ## helper function to process absolute abundance data
  process_absolute <- function(absolute, llod_inds, abs_inds, llod = 0, div_num = 100) {
    ## subset with absolute abundance values
    sub_absolute <- absolute[, abs_inds]/div_num
    sub_llod <- absolute[, llod_inds]/div_num
    ## clean up with llod
    clean_absolute <- mapply(function(x, y) ifelse(x < y, llod, x), sub_absolute, sub_llod)
    return(clean_absolute)
  }
  ## if full_data is a matrix, make it a data.frame
  if (is.matrix(full_data)) {
    tmp <- as.data.frame(full_data)
    full_data <- tmp
  }

  ## subset the absolute abundance data; set counts at lower limit of detection to llod
  sub_absolute <- full_data[, abs_inds]
  clean_absolute <- process_absolute(sub_absolute, grepl(regex_thr, names(sub_absolute), fixed = TRUE), grepl(regex_abs, names(sub_absolute), fixed = TRUE), llod, div_num)
  tmp_absolute <- data.frame(clean_absolute)

  ## subset the relative abundance data
  subset_relative <- full_data[, rel_inds]
  ## re-order so that the bugs with absolute abundance are first
  all_relative <- cbind(subset_relative[, abs_plus_rel_inds], subset_relative[, -c(abs_plus_rel_inds)])
  ## check whether or not any have a nonzero count
  nonzero_counts <- apply(all_relative > 0, 2, sum) > 0
  tmp_relative <- data.frame(all_relative[, nonzero_counts])
  names(tmp_relative) <- names(full_data)[rel_inds]

  ## remove any rows that have less than m_min reads
  m <- apply(tmp_relative, 1, sum)
  absolute_m <- tmp_absolute[m >= m_min, ]
  relative_m <- tmp_relative[m >= m_min, ]
  ## remove any rows that have missing absolute abundances
  na_absolute <- apply(absolute_m, 1, function(x) any(is.na(x)))
  absolute <- absolute_m[!na_absolute, ]
  relative <- relative_m[!na_absolute, ]
  return(list(absolute = absolute, relative = relative))
}
