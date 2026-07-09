biclust <- function(
  mat,
  k,
  use_test = TRUE,
  n_reps = 20,
  tol = 0.01,
  verbose = FALSE,
  ...
) {

  ## 1. Build consensus
  if (verbose){
    message("Starting consensus NMF.")
  }

  result <- .build_consensus_nmf(
    mat$train,
    k,
    reps = n_reps,
    tol = tol,
    ...
  )


  ## 2. Project test data if available
  if (use_test) {

    if (verbose){
      message("Projecting on held-out test data.")
    }

    ## only if mat actually contains test data
    result <- .project_on_test(
      result,
      test_data = mat$test
    )
  }

  ## 3. Extract factors

  result <- .get_factors(result)
  if (verbose){
    message("Consensus reached.")
  }
  ## 4. Extract candidate biclusters

  result <- .extract_core(result)
  if (verbose){
    message("Candidate bicluster found")
  }
  result
}