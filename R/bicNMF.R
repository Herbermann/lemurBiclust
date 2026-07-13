#' @export
bicNMF <- function(
  mat,
  k,
  n_reps = 20,
  tol = 0.01,
  verbose = FALSE,
  threads = 1,
  threads_RcppML = 1,
  L1 = c(0.0, 0.0),
  zero_tol = 1e-6,
  backend = c("serial", "parallel"),
  gene_support_thrs = 0.9,
  cell_support_thrs = 0.9
) {

  if (verbose){
    message("Starting consensus NMF.")
  }

  result <- .build_consensus_nmf(
    mat,
    k,
    reps = n_reps,
    tol = tol,
    verbose = verbose,
    threads = threads,
    threads_RcppML = threads_RcppML,
    L1 = L1,
    zero_tol = zero_tol,
    backend = backend
  )


  result <- .get_factors(result, gene_support_thrs, cell_support_thrs)
  if (verbose){
    message("Consensus reached.")
  }

  result <- .extract_core(result)
  if (verbose){
    message("Candidate bicluster found")
  }
  result
}



