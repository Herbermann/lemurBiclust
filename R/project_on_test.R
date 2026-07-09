.project_on_test <- function(
  res,
  test_data
) {
  
  consensus_obj <- res$consensus

  n_consensus <- length(consensus_obj$models)
  zero_tol <- consensus_obj$metadata$zero_tol

  # Project each fitted model onto the test data, then align factor order
  list_H_test <- lapply(seq_len(n_consensus), function(i) {

    current_model <- consensus_obj$models[[i]]
    H_new <- RcppML::predict(current_model, test_data)@h

    # normalize rows safely
    rs <- rowSums(H_new, na.rm = TRUE)
    rs[rs == 0] <- 1
    H_new <- H_new / rs

    # align test factors back to the reference order
    map <- consensus_obj$alignments[[i]]$map

    H_aligned <- matrix(
      0,
      nrow = length(map),
      ncol = ncol(H_new),
      dimnames = list(
        paste0("factor_", seq_len(length(map))),
        colnames(H_new)
      )
    )

    for (f in seq_len(length(map))) {
      mf <- map[f]
      H_aligned[f, ] <- H_new[mf, ]
    }

    H_aligned
  })

  k <- nrow(list_H_test[[1]])
  n_samples <- ncol(list_H_test[[1]])
  sample_names <- colnames(list_H_test[[1]])
  if (is.null(sample_names)) {
    sample_names <- paste0("sample_", seq_len(n_samples))
  }

  cell_test_votes <- matrix(
    0,
    nrow = n_samples,
    ncol = k,
    dimnames = list(sample_names, paste0("factor_", seq_len(k)))
  )

  cell_test_loading_sum <- matrix(
    0,
    nrow = n_samples,
    ncol = k,
    dimnames = list(sample_names, paste0("factor_", seq_len(k)))
  )

  for (r in seq_along(list_H_test)) {

    H <- abs(list_H_test[[r]])

    for (f in seq_len(k)) {

      cell_idx <- which(H[f, ] > zero_tol)

      if (length(cell_idx)) {
        cell_test_votes[cell_idx, f] <- cell_test_votes[cell_idx, f] + 1
        cell_test_loading_sum[cell_idx, f] <-
          cell_test_loading_sum[cell_idx, f] + H[f, cell_idx]
      }
    }
  }

  cell_test_support <- cell_test_votes / n_consensus
  cell_test_loading_mean <- cell_test_loading_sum / pmax(cell_test_votes, 1)

  consensus_obj$cell_test_support <- cell_test_support
  consensus_obj$cell_test_votes <- cell_test_votes
  consensus_obj$cell_test_loading_mean <- cell_test_loading_mean
  consensus_obj$test_names <- sample_names
  consensus_obj$h_test <- list_H_test

  res$consensus <- consensus_obj

  return(res)
}






project_on_test <- function(consensus_obj, test_data) {

  n_consensus <- length(consensus_obj$models)
  zero_tol <- 1e-6

  # Project each fitted model onto the test data, then align factor order
  list_H_test <- lapply(seq_len(n_consensus), function(i) {

    current_model <- consensus_obj$models[[i]]
    H_new <- RcppML::predict(current_model, test_data)@h

    # normalize rows safely
    rs <- rowSums(H_new, na.rm = TRUE)
    rs[rs == 0] <- 1
    H_new <- H_new / rs

    # align test factors back to the reference order
    map <- consensus_obj$alignments[[i]]$map

    H_aligned <- matrix(
      0,
      nrow = length(map),
      ncol = ncol(H_new),
      dimnames = list(
        paste0("factor_", seq_len(length(map))),
        colnames(H_new)
      )
    )

    for (f in seq_len(length(map))) {
      mf <- map[f]
      H_aligned[f, ] <- H_new[mf, ]
    }

    H_aligned
  })

  k <- nrow(list_H_test[[1]])
  n_samples <- ncol(list_H_test[[1]])
  sample_names <- colnames(list_H_test[[1]])
  if (is.null(sample_names)) {
    sample_names <- paste0("sample_", seq_len(n_samples))
  }

  cell_test_votes <- matrix(
    0,
    nrow = n_samples,
    ncol = k,
    dimnames = list(sample_names, paste0("factor_", seq_len(k)))
  )

  cell_test_loading_sum <- matrix(
    0,
    nrow = n_samples,
    ncol = k,
    dimnames = list(sample_names, paste0("factor_", seq_len(k)))
  )

  for (r in seq_along(list_H_test)) {

    H <- abs(list_H_test[[r]])

    for (f in seq_len(k)) {

      cell_idx <- which(H[f, ] > zero_tol)

      if (length(cell_idx)) {
        cell_test_votes[cell_idx, f] <- cell_test_votes[cell_idx, f] + 1
        cell_test_loading_sum[cell_idx, f] <-
          cell_test_loading_sum[cell_idx, f] + H[f, cell_idx]
      }
    }
  }

  cell_test_support <- cell_test_votes / n_consensus
  cell_test_loading_mean <- cell_test_loading_sum / pmax(cell_test_votes, 1)

  consensus_obj$cell_test_support <- cell_test_support
  consensus_obj$cell_test_votes <- cell_test_votes
  consensus_obj$cell_test_loading_mean <- cell_test_loading_mean
  consensus_obj$test_names <- sample_names
  consensus_obj$h_test <- list_H_test

  consensus_obj
}