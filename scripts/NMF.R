devtools::load_all("/Users/herbermann/code/lemur_biclusters/lemur")

library("lemur")
library("tidyverse")
library("SingleCellExperiment")
library("patchwork")
library("RcppML")
package.version("RcppML")

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


#############################



estimate_factors <- function(){
    
  ks <- c(2, 3, 4, 5, 6, 8, 10, 12, 14)
  n_restarts <- 10

  results <- lapply(ks, function(k) {
    # multiple restarts, keep best
    fits <- lapply(seq_len(n_restarts), function(i) {
      nmf(de_pos, k = k, seed = i, tol = 1e-4)
    })
    errors <- sapply(fits, function(f) {
      sum((de_pos - f$w %*% diag(f$d) %*% f$h)^2)
    })
    best <- fits[[which.min(errors)]]
    list(k = k, error = min(errors), fit = best)
  })

  # plot reconstruction error vs k
  plot(ks, sapply(results, `[[`, "error"), type = "b",
      xlab = "k", ylab = "reconstruction error")

}


estimate_stability <- function(){
# cosine similarity between W matrices across restarts
# for a given k, high mean cosine similarity = stable factors
  n_restarts <- 10
  ks <- c(2, 3, 4, 5, 6, 8, 10, 12, 14)

  stability <- sapply(ks, function(k) {
  fits <- lapply(seq_len(n_restarts), function(i) {
    nmf(de_pos, k = k, seed = i, tol = 1e-4)
  })
  ws <- lapply(fits, function(m) {
    w <- m@w
    live <- colSums(w) > 0
    w[, live, drop = FALSE]
  })  # pairwise cosine similarity of factor columns across pairs of fits
  sims <- combn(length(ws), 2, function(idx) {
    W1 <- ws[[idx[1]]]
    W2 <- ws[[idx[2]]]
    # match factors by max cosine similarity (greedy)
    W1 <- scale(W1, center = FALSE, scale = sqrt(colSums(W1^2)))
    W2 <- scale(W2, center = FALSE, scale = sqrt(colSums(W2^2)))
    sim_mat <- t(W1) %*% W2
    mean(apply(sim_mat, 1, max))
  }, simplify = TRUE)
  mean(sims)
})

plot(ks, stability, type = "b",
     xlab = "k", ylab = "mean factor stability")
}



k <- 8
L1 <- c(0.5, 0.1)

fit_pos <- nmf(de_pos, k = k, L1 = L1, 
               seed = 1, tol = 1e-5)

fit_neg <- nmf(de_neg, k = k, L1 = L1, 
               seed = 1, tol = 1e-5)

W <- fit_pos$w
W_norm <- scale(W, center = FALSE, scale = sqrt(colSums(W^2)))
round(t(W_norm) %*% W_norm, 3)



plot_distributions <- function(){
  # distribution of cell scores for each factor
  par(mfrow = c(2,4))
  for(i in 1:k) {
    hist(fit_pos$h[i,], breaks = 50, 
        main = paste("pos factor", i),
        xlab = "cell score")
    hist(fit_neg$h[i,], breaks = 50,
        main = paste("neg factor", i), 
        xlab = "cell score")
  }

  par(mfrow = c(2,4))
  for(i in 1:k) {
    hist(fit_pos$w[,i], breaks = 50, 
        main = paste("pos factor", i),
        xlab = "gene score")
    hist(fit_neg$w[,i], breaks = 50,
        main = paste("neg factor", i), 
        xlab = "gene score")
  } 
}

generate_null_dist <- function(matrix, k, L1 = c(0, 0), n_perm = 20) {
  
  nulls <- lapply(seq_len(n_perm), function(i) {
    mat_perm <- t(apply(matrix, 1, sample))
    fit_perm <- nmf(mat_perm, k = k, L1 = L1, seed = i, tol = 1e-4)
    list(h = fit_perm@h, w = fit_perm@w)
  })

  null_h <- lapply(nulls, `[[`, "h")
  null_w <- lapply(nulls, `[[`, "w")

  null_threshold_h <- sapply(seq_len(k), function(i) {
    quantile(unlist(lapply(null_h, function(h) h[i, ])), 0.95)
  })

  null_threshold_w <- sapply(seq_len(k), function(i) {
    quantile(unlist(lapply(null_w, function(w) w[, i])), 0.95)
  })

  list(h = null_threshold_h, w = null_threshold_w)
}

null_pos_ <- generate_null_dist(de_pos, k = k, L1 = L1)
null_neg_ <- generate_null_dist(de_neg, k = k, L1 = L1)

null_pos   <- null_pos_$h
null_neg   <- null_neg_$h
null_pos_W <- null_pos_$w
null_neg_W <- null_neg_$w


cell_regions <- tibble(
  umap1 = umap[, 1],
  umap2 = umap[, 2]
)
for (i in seq_len(nrow(fit_pos@h))) {
  cell_regions[[paste0("program_", LETTERS[i], "p")]] <- fit_pos@h[i, ] > null_pos[i]
}
for (i in seq_len(nrow(fit_neg@h))) {
  cell_regions[[paste0("program_", LETTERS[i], "n")]] <- fit_neg@h[i, ] > null_neg[i]
}

cell_regions |>
  pivot_longer(starts_with("program_"),
               names_to = "program", values_to = "score") |>
  ggplot(aes(x = umap1, y = umap2, color = score)) +
    geom_point(size = 0.3) +
    scale_color_gradient2(low = "#FFD800", high= "#0056B9") +
    facet_wrap(~ program) +
    coord_fixed()



jaccard_with_indices <- function(A, B){
  num <- length(intersect(A, B))
  denom <- length(union(A, B))
  num/denom
}

jaccard_with_logic <- function(A, B){
  sum(A & B) / sum(A | B)
}

get_jaccard_matrix_of_regions <- function(cell_regions){

  programs <- cell_regions[, grepl("^program_", names(cell_regions))]
  jaccard_matrix <- outer(
    names(programs),
    names(programs),
    Vectorize(function(a, b) jaccard_with_logic(programs[[a]], programs[[b]]))
  )

  rownames(jaccard_matrix) <- names(programs)
  colnames(jaccard_matrix) <- names(programs)

  jaccard_matrix
}



genespace_stability <- function(W1, W2) {
  # principal angles between column spaces of W1 and W2
  # high similarity = same subspace even if individual factors differ
  Q1 <- qr.Q(qr(W1))
  Q2 <- qr.Q(qr(W2))
  decomp <- svd(t(Q1) %*% Q2)
  decomp@d
}
