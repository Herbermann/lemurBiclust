#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
library(dplyr)
library(data.table)
library(SingleCellExperiment)
library(Matrix)
library(lemur)
library(tidyverse)
library(Polychrome)
library(scico)
library(patchwork)
library("progressr")
set.seed(42)
#
#
#
#
#
folder_data <- "/Users/herbermann/code/data/coculture/"
file_name <- "coculture_effect_in_AML_5000"
wd <- "/Users/herbermann/code/lemur_biclusters"
cache_dir <- file.path(wd, "rds")
#
#
#
lemur_fit_file <- file.path(cache_dir, "lemur_essi_fit.rds")

FLAG_fit_lemur <- TRUE

if (FLAG_fit_lemur == TRUE){
  if (file.exists(lemur_fit_file)) {
    fit <- readRDS(lemur_fit_file)
  } else {
    
    sce_all <- readRDS(paste0(folder_data, "sce_object_", file_name, ".RDS"))
    
    n_latents <- 20
    fit <- lemur(sce_all, ~coculture, n_embedding = n_latents, test_fraction = 0.5)
    fit <- align_by_grouping(fit, grouping = colData(fit)[, "clusters_renamed"])
    fit <- test_de(fit, contrast = cond(coculture = "Coculture") - cond(coculture = "BM_only"))
    
    saveRDS(fit, lemur_fit_file)
  }
}
#
#
#
#
#
#
#
#
#
library(Matrix)
library("RcppML")
package.version("RcppML")

de_matrix <- as.matrix(assay(fit$training_data, "DE"))

de_pos <- pmax(de_matrix, 0.)   # genes x cells, upregulation only
de_neg <- pmax(-de_matrix, 0.)  # genes x cells, downregulation only

quant_pos <- quantile(de_pos, c(0.5, 0.75, 0.9, 0.95))
quant_neg <- quantile(de_neg, c(0.5, 0.75, 0.9, 0.95))

de_pos <- tanh(2/quant_pos[4] * de_pos)
de_neg <- tanh(2/quant_neg[4] * de_neg)

de_pos <- Matrix(de_pos, sparse=TRUE)
de_neg <- Matrix(de_neg, sparse=TRUE)

pos_orig <- rownames(de_pos)
neg_orig <- rownames(de_neg)

rownames(de_pos) <- paste0("up_", make.unique(pos_orig))
rownames(de_neg) <- paste0("down_", make.unique(neg_orig))

de <- rbind(de_pos, de_neg)

#
#
#
#
#
#

genespace_stability <- function(fit1, fit2) {
  # weight W columns by their factor importance
  W1 <- fit1@w %*% diag(fit1@d)
  W2 <- fit2@w %*% diag(fit2@d)
  
  # remove dead factors
  W1 <- W1[, colSums(W1) > 0, drop = FALSE]
  W2 <- W2[, colSums(W2) > 0, drop = FALSE]
  
  # now QR on the weighted W
  Q1 <- qr.Q(qr(W1))
  Q2 <- qr.Q(qr(W2))
  
  decomp <- svd(t(Q1) %*% Q2)
  theta2 <- acos(decomp@d)^2
  sqrt(theta2)/(length(theta2)* 0.5 * 3.1415)
}

find_factors <- function(
    matrix,
    ks,
    n_restarts = 2,
    test_fraction = 0.,
    L1 = c(0., 0.),
    tol = 1e-4,
    show_plot = FALSE,
    stability_thrsh = 0.9 
){
  
  p <- progressr::progressor(steps = length(ks) * n_restarts)

  results <- lapply(ks, function(k) {
    
    # single loop over restarts, compute both error and stability from same fits
    if (test_fraction == 0.){
      fits <- lapply(seq_len(n_restarts), function(i) {
      p(message = sprintf("k=%d, restart %d/%d", k, i, n_restarts))
      nmf(matrix, k = k, seed = i, L1 = L1, tol = tol)
      })
    }
    else {
      fits <- lapply(seq_len(n_restarts), function(i) {
      p(message = sprintf("k=%d, restart %d/%d", k, i, n_restarts))
      nmf(matrix, k = k, seed = i, L1 = L1, tol = tol, test_fraction = test_fraction, mask="zeros")
      })
    }
    
    # --- reconstruction error ---
    loss_slot <- if (test_fraction == 0) "loss" else "test_loss"
    errors <- sapply(fits, function(f) f@misc[[loss_slot]])
    best_idx <- which.min(errors)
    
    # --- stability ---
    ws <- lapply(fits, function(m) {
      w <- m@w
      w[, colSums(w) > 0, drop = FALSE]
    })

    sims <- combn(length(fits), 2, function(idx) {
      if (ncol(fits[[idx[1]]]@w) == 0 || ncol(fits[[idx[2]]]@w) == 0) return(NA_real_)
      d <- genespace_stability(fits[[idx[1]]], fits[[idx[2]]])
      mean(d)
    }, simplify = TRUE)

    out <- list(
      k        = k,
      error    = errors[best_idx],
      stability= mean(sims, na.rm = TRUE)
      #fit      = fits[[best_idx]]
    )

    out
  })

  error_df <- data.frame(
    k     = ks,
    error = sapply(results, `[[`, "error")
  )
  stab_df <- data.frame(
    k         = ks,
    stability = sapply(results, `[[`, "stability")
  )

  p_error <- ggplot(error_df, aes(x = k, y = error)) +
    geom_line() + geom_point() +
    labs(x = "k", y = "loss or reconstr. error") +
    theme_classic()

  p_stab <- ggplot(stab_df, aes(x = k, y = stability)) +
    geom_line() + geom_point() +
    labs(x = "k", y = "factor stability (PA spans gene space)") +
    theme_classic()

  pl <- (p_error | p_stab)

  return(list(
    results   = results,
    error_df  = error_df,
    stab_df   = stab_df,
    plot      = pl
  ))
}


filename <- file.path(cache_dir, "find_k_coarse.RDS")
if (file.exists(filename)) {
  find_k_coarse <- readRDS(filename)
} else {
    find_k_coarse <- find_factors(
      de,
      ks = c(8, 10, 12, 14, 16, 18),
      n_restarts=3,
      L1 = c(0.0 , 0.0),
      tol = 1e-2
    )
    saveRDS(find_k_coarse, filename)
}


filename <- file.path(cache_dir, "find_k_fine.RDS")
if (file.exists(filename)) {
  find_k_fine <- readRDS(filename)
} else {
    find_k_fine <- find_factors(
      de,
      ks = c(14, 15, 16),
      n_restarts=3,
      L1 = c(0.0 , 0.0),
      tol = 1e-4,
      test_fraction = 0.1
    )
    saveRDS(find_k_fine, filename)
}
#
#
#
#
#
fit_nmf <- function(matrix,
              k,
              L1 = c(0, 0),
              n_perm = 20,
              thr_h = 0.95,
              thr_w_pos = 0.95,
              thr_w_neg = 0.95
              )
  {

  master_result <- nmf()
  
  n_genes <- nrow(matrix) / 2

  nulls <- lapply(seq_len(n_perm), function(i) {
    mat_perm <- t(apply(matrix, 1, sample))
    fit_perm <- nmf(mat_perm, k = k, L1 = L1, seed = i, tol = 1e-4)
    list(h = fit_perm@h, w = fit_perm@w)
  })

  null_h <- lapply(nulls, `[[`, "h")
  null_w <- lapply(nulls, `[[`, "w")

  null_w_pos <- lapply(null_w, function(w) w[1:n_genes, , drop = FALSE])
  null_w_neg <- lapply(null_w, function(w) w[(n_genes + 1):nrow(w), , drop = FALSE])

  null_threshold_h <- sapply(seq_len(k), function(i) {
    quantile(unlist(lapply(null_h, function(h) h[i, ])), thr_h)
  })

  null_threshold_w_pos <- sapply(seq_len(k), function(i) {
    quantile(unlist(lapply(null_w_pos, function(w) w[, i])), thr_w_pos)
  })

  null_threshold_w_neg <- sapply(seq_len(k), function(i) {
    quantile(unlist(lapply(null_w_neg, function(w) w[, i])), thr_w_neg)
  })

  list(h = null_threshold_h, w_pos = null_threshold_w_pos, w_neg = null_threshold_w_neg)
}
#
#
#
