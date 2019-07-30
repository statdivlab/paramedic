#' Process a dataset, returning qPCR and br16S matrices
#'
#' Return matrices of observed qPCR and br16S abundances, given a single dataset and a set of indices.
#'
#' @param full_data the full dataset, containing both qPCR and br16S abundances
#' @param br_inds the indices marking the br16S measurements in \code{full_data}
#' @param qpcr_inds the indices marking the qPCR measurements in \code{full_data}
#' @param pcr_plus_br_inds the indices that have both br16S and qPCR measurements
#' @param regex_thr the regular expression to identify qPCR threshold for each taxon-specific qPCR (e.g., "_t" appended to each threshold)
#' @param regex_cps the regular expression to identify qPCR measurement for each taxon-specific qPCR (e.g., "_cps" appended to each measurement of copies per swab)
#' @param llod the lower limit of detection for qPCR
#' @param m_min the minimum number of reads to consider
#' @param div_num the number to divide qPCR by (to create valid R integers for use in Stan)
#'
#' @return a list containing two matrices: one with the observed qPCR, and the other with the observed br16S.
#'
#' @details
#'
#' @examples
#' ## load the package, read in example data
#' library("paramedic")
#' data(simple_example_data)
#'
#' ## process the example data
#' q <- 3
#' q_obs <- 2
#' processed_data <- process_data(full_data = simple_example_data, br_inds = 1:q,
#' qpcr_inds = (q + 1):(q + 1 + q_obs),
#' pcr_plus_br_inds = 1:q_obs,
#' regex_thr = NA, regex_cps = "", llod = 0,
#' m_min = 1000, div_num = 1)
#'
#' @export
process_data <- function(full_data, br_inds, qpcr_inds, pcr_plus_br_inds, regex_thr = "_t", regex_cps = "_cps", llod = 0, m_min = 1000, div_num = 100) {
  ## helper function to process qPCR
  process_qpcr <- function(qpcr, llod_inds, qpcr_inds, llod = 0, div_num = 100) {
    ## subset with qpcr values
    sub_qpcr <- qpcr[, qpcr_inds]/div_num
    sub_llod <- qpcr[, llod_inds]/div_num
    ## clean up with llod
    clean_qpcr <- mapply(function(x, y) ifelse(x < y, llod, x), sub_qpcr, sub_llod)
    return(clean_qpcr)
  }
  ## if full_data is a matrix, make it a data.frame
  if (is.matrix(full_data)) {
    tmp <- as.data.frame(full_data)
    full_data <- tmp
  }

  ## subset the qpcr data; set counts at lower limit of detection to llod
  sub_qpcr <- full_data[, qpcr_inds]
  clean_qpcr <- process_qpcr(sub_qpcr, grepl(regex_thr, names(sub_qpcr), fixed = TRUE), grepl(regex_cps, names(sub_qpcr), fixed = TRUE), llod, div_num)
  tmp_qpcr <- data.frame(clean_qpcr)

  ## subset the br16s data
  subset_br <- full_data[, br_inds]
  ## re-order so that the bugs with qpcr are first
  all_br <- cbind(subset_br[, pcr_plus_br_inds], subset_br[, -c(pcr_plus_br_inds)])
  ## check whether or not any have a nonzero count
  nonzero_counts <- apply(all_br > 0, 2, sum) > 0
  tmp_br <- data.frame(all_br[, nonzero_counts])
  names(tmp_br) <- names(full_data)[br_inds]

  ## remove any rows that have less than m_min reads
  m <- apply(tmp_br, 1, sum)
  qpcr_m <- tmp_qpcr[m >= m_min, ]
  br_m <- tmp_br[m >= m_min, ]
  ## remove any rows that have missing qpcr
  na_qpcr <- apply(qpcr_m, 1, function(x) any(is.na(x)))
  qpcr <- qpcr_m[!na_qpcr, ]
  br <- br_m[!na_qpcr, ]
  return(list(qpcr = qpcr, br = br))
}
