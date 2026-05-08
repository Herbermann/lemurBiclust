.project_biclusters <- function(
    best_fit,
    test_mat,
    thresholds = NULL,
    train_H = best_fit$H,
    return_thresholded = TRUE
) {

  # ------------------------------------------------------------
  # Predict H on held-out data
  # ------------------------------------------------------------

  H_test <- .predict_nmf(best_fit, test_mat)

  rownames(H_test) <- rownames(train_H)
  colnames(H_test) <- colnames(test_mat)

  k <- nrow(H_test)

  bic_names <- rownames(H_test)
  if (is.null(bic_names)) {
    bic_names <- paste0("Bic", seq_len(k))
    rownames(H_test) <- bic_names
  }

  out <- list(
    H_test = H_test
  )

  # ------------------------------------------------------------
  # Optional thresholding / gatekeeping
  # ------------------------------------------------------------

  if (!is.null(thresholds) && return_thresholded) {

    stopifnot(length(thresholds) == k)

    thresholded <- lapply(seq_len(k), function(i) {

      vals <- sort(H_test[i, ], decreasing = TRUE)

      .fast_extract(vals, thresholds[i])
    })

    names(thresholded) <- bic_names

    out$thresholds <- thresholds
    out$thresholded <- thresholded
  }

  out
}