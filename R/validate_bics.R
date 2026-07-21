.validate_overlap_bic <- function(
  core_object,
  embedding_t,
  j
) {


  cells_A <- core_object$train[[j]]
  cells_B <- core_object$test[[j]]

  cells <- c(cells_A, cells_B)
  X <- embedding_t[cells, , drop = FALSE]
  group <- factor(ifelse(cells %in% cells_A, "A", "B"))

  n <- min(length(cells_A), length(cells_B))
  k <- min(floor(sqrt(n)), 50, n - 1)

  nn <- FNN::get.knn(X, k = k)$nn.index
  nn_group <- matrix(group[nn], nrow = nrow(nn))
  mix_per_cell <- rowMeans(nn_group != group)

  mix_score <- mean(mix_per_cell)

  nA <- length(cells_A)
  nB <- length(cells_B)

  pA <- nA / (nA + nB)
  expected <- 2 * pA * (1 - pA)

  mix_score_norm <- mix_score / expected

  pA_local <- rowMeans(nn_group == "A")
  pB_local <- 1 - pA_local

  lisi_per_cell <- 1 / (pA_local^2 + pB_local^2)
  mean_lisi <- mean(lisi_per_cell)

  return(c("kNN_mixing" = mix_score_norm, "LISI" = mean_lisi))
}


validate_overlap <- function(
  core_object,
  embedding
) {

  n_bic <- length(core_object$train)
  results <- data_frame(
    bicluster = seq(n_bic),
    kNN_mixing = numeric(n_bic),
    LISI = numeric(n_bic)
  )

  embedding_t <- t(embedding)

  for (i in seq(n_bic)){
    temp <- .validate_overlap_bic(core_object, embedding_t, i)
    results$kNN_mixing[i] <- temp[1]
    results$LISI[i] <- temp[2]
  }

  results
}
