devtools::load_all("/Users/herbermann/code/lemur_biclusters/lemur")
library("RcppML")
#######

library("lemur")
library("tidyverse")
library("SingleCellExperiment")
library("patchwork")
library("RcppML")
library("progressr")


set.seed(42)

data("glioblastoma_example_data", package = "lemur")
glioblastoma_example_data

fit <- lemur(glioblastoma_example_data, design = ~ patient_id + condition, 
             n_embedding = 15, test_fraction = 0.5)

fit <- align_harmony(fit)

umap <- uwot::umap(t(fit$embedding))

fit <- compute_contrasts(fit, contrast = cond(condition = "panobinostat") - cond(condition = "ctrl"))

sel_gene <- "ENSG00000169429" # is CXCL8

p1 <- tibble(umap = umap) |>
  mutate(de = assay(fit, "DE")[sel_gene,]) |>
  arrange(abs(de)) |>
  ggplot(aes(x = umap[,1], y = umap[,2])) +
    geom_point(aes(color = de), size = 0.5) +
    scale_color_gradient2(low = "#FFD800", high= "#0056B9") +
    coord_fixed()


######################
# LEGACY NEIGHBORHOODS

neighborhoods <- find_de_neighborhoods(fit, group_by = vars(patient_id, condition))

as_tibble(neighborhoods) |>
  left_join(as_tibble(rowData(fit)[,1:2]), by = c("name" = "gene_id")) |>
  relocate(symbol, .before = "name") |>
  arrange(pval) |>
  head(5)


##################################
# ------ PREPARE MATRIX ---------#
##################################

library(Matrix)

de_matrix <- as.matrix(assay(fit, "DE"))

de_pos <- pmax(de_matrix, 0.2)   # genes x cells, upregulation only
de_neg <- pmax(-de_matrix, 0.2)  # genes x cells, downregulation only

quant_pos <- quantile(de_pos, c(0.5, 0.75, 0.9, 0.95))
quant_neg <- quantile(de_neg, c(0.5, 0.75, 0.9, 0.95))

de_pos <- tanh(2/quant_pos[4] * de_pos)
de_neg <- tanh(2/quant_neg[4] * de_neg)

de_pos <- Matrix(de_pos, sparse=TRUE)
de_neg <- Matrix(de_neg, sparse=TRUE)



##################################



##################################



##################################

#genespace_stability <- function(W1, W2) {
#  # principal angles between column spaces of W1 and W2
#  # high similarity = same subspace even if individual factors differ
#  Q1 <- qr.Q(qr(W1))
#  Q2 <- qr.Q(qr(W2))
#  decomp <- svd(t(Q1) %*% Q2)
#  #mean(decomp@d)
#  theta2 <- acos(decomp@d)^2
#  sqrt(theta2)/(length(theta2)*0.5 * 3.1415)
#}

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
      stability = mean(sims, na.rm = TRUE)
      #fit      = fits[[best_idx]]
    )
    out
  })


  if (show_plot == TRUE){

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
      labs(x = "k", y = "reconstruction error") +
      theme_classic()

    p_stab <- ggplot(stab_df, aes(x = k, y = stability)) +
      geom_line() + geom_point() +
      labs(x = "k", y = "mean factor stability") +
      theme_classic()

    print(p_error + p_stab)

  }

  results

}



with_progress({
  test <- find_factors(de_pos, c(8, 10, 12, 14, 16), n_restarts=3, test_fraction = 0.1, show_plot = TRUE, L1 = c(0.0 , 0.0) )
})


l1_sweep <- c(0, 0.01, 0.05, 0.1, 0.2)
fits_l1 <- lapply(l1_sweep, function(l1) {
  nmf(de_pos, k = 12, seed = 1:5, L1 = c(l1, 0), tol = 1e-4)
})

# track mean off-diagonal cosine overlap as L1 increases
overlaps <- sapply(fits_l1, function(f) {
  C <- cosine(f@w, f@w)
  diag(C) <- NA
  mean(C, na.rm = TRUE)
})

data.frame(L1 = l1_sweep, overlap = overlaps)

mean_overlap <- function(fit, top_n = 50) {
  W <- fit@w
  # for each factor, get indices of top n genes by loading
  top_genes <- apply(W, 2, function(w) order(w, decreasing = TRUE)[seq_len(top_n)])
  
  # pairwise Jaccard similarity of top gene sets
  pairs <- combn(ncol(W), 2, function(idx) {
    g1 <- top_genes[, idx[1]]
    g2 <- top_genes[, idx[2]]
    length(intersect(g1, g2)) / length(union(g1, g2))
  }, simplify = TRUE)
  
  mean(pairs)
}





imputation_cv <- function(matrix,
                          ks,
                          n_restarts = 10,
                          mask_fraction = 0.3,
                          L1 = c(0, 0),
                          tol = 1e-4) {
  
  p <- progressr::progressor(steps = length(ks) * n_restarts)
  
  results <- lapply(ks, function(k) {
    
    losses <- sapply(seq_len(n_restarts), function(i) {
      p()
      # mask only from non-zero entries
      nonzero_idx <- which(matrix > 0, arr.ind = TRUE)
      selected    <- nonzero_idx[sample(nrow(nonzero_idx), 
                                        floor(mask_fraction * length(nonzero_idx))), ]
      mask        <- Matrix::sparseMatrix(
        i = selected[, 1], j = selected[, 2],
        dims = dim(matrix), x = TRUE
      )
      mask <- as(mask, "dgCMatrix")

      # fit with mask — masked entries excluded from ALS
      fit <- nmf(matrix, k = k, seed = i, L1 = L1, tol = tol, mask = mask)
      
      # evaluate loss only on masked entries
      reconstructed <- fit@w %*% diag(fit@d) %*% fit@h
      masked_idx <- which(mask > 0)
      mean((matrix[masked_idx] - reconstructed[masked_idx])^2)
    })
    
    list(k = k, mean_loss = mean(losses), sd_loss = sd(losses))
  })
  
  df <- data.frame(
    k         = sapply(results, `[[`, "k"),
    mean_loss = sapply(results, `[[`, "mean_loss"),
    sd_loss   = sapply(results, `[[`, "sd_loss")
  )
  
  p_out <- ggplot(df, aes(x = k, y = mean_loss)) +
    geom_line() + geom_point() +
    geom_ribbon(aes(ymin = mean_loss - sd_loss,
                    ymax = mean_loss + sd_loss), alpha = 0.2) +
    labs(x = "k", y = "mean imputation MSE") +
    theme_classic()
  
  list(plot = p_out, data = df)
}

with_progress({
  test <- imputation_cv(de_pos, c(10, 11, 12), n_restarts=10, mask_fraction = 0.3) 
})
test$plot


biclusters <- nmf(de_pos, k = 12, seed = 1:10)
W <- biclusters@w
H <- biclusters@h

rownames(W) <- rownames(assay(fit, "DE"))
colnames(H) <- colnames(assay(fit, "DE"))

hist(W[,12])
hist(H[,12])


###############
generate_null_dist <- function(matrix, k, L1 = c(0, 0), n_perm = 20, thr_h=0.95, thr_w=0.95) {
  
  nulls <- lapply(seq_len(n_perm), function(i) {
    mat_perm <- t(apply(matrix, 1, sample))
    fit_perm <- nmf(mat_perm, k = k, L1 = L1, seed = i, tol = 1e-4)
    list(h = fit_perm@h, w = fit_perm@w)
  })

  null_h <- lapply(nulls, `[[`, "h")
  null_w <- lapply(nulls, `[[`, "w")

  null_threshold_h <- sapply(seq_len(k), function(i) {
    quantile(unlist(lapply(null_h, function(h) h[i, ])), thr_h)
  })

  null_threshold_w <- sapply(seq_len(k), function(i) {
    quantile(unlist(lapply(null_w, function(w) w[, i])), thr_w)
  })

  list(h = null_threshold_h, w = null_threshold_w)
}

null_pos_ <- generate_null_dist(de_pos, k = 12, n_perm = 10)

null_pos   <- null_pos_$h
null_pos_W <- null_pos_$w


#################
# Plottings.

k <- 11
l <- 2

de_expression <- assay(fit, "DE")

genes <- sort(W[,k], decreasing = TRUE)
genes[l]
sel_gene <- names(genes[l])

p1 <- tibble(umap = umap) |>
  mutate(de = assay(fit, "DE")[sel_gene,]) |>
  #arrange(abs(de)) |>
  ggplot(aes(x = umap[,1], y = umap[,2])) +
    geom_point(aes(color = de), size = 0.5) +
    scale_color_gradient2(low = "#FFD800", high= "#0056B9") +
    coord_fixed()


null_thresh_h <- null_pos[k]
null_thresh_w <- null_pos_W[k]

member_cells <- H[k,]>null_thresh_h
member_genes <- W[,k]>null_thresh_w


p2 <- tibble(umap = umap) |>
  mutate(member_cells = member_cells) |>
  ggplot(aes(x = umap[,1], y = umap[,2])) +
    geom_point(aes(color = member_cells), size = 0.5) +
    coord_fixed()

p1 | p2


out_nei <- de_expression[sel_gene, which(member_cells == FALSE)]
in_nei <- de_expression[sel_gene, which(member_cells == TRUE)]

t.test(in_nei, out_nei)

