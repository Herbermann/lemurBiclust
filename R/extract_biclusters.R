.predict_nmf <- function(fit, mat) {
  # We have to correct manually for predict returning unnormaised loadings...
  pred <- RcppML::predict(fit$fit, mat)
  sweep(pred@h, 1, fit$D, "/")
}

.null_pvalue <- function(x, null_vec, probs) {
  1 - approx(x = null_vec, y = probs, xout = x, rule = 2)$y
}

extract_biclusters_legacy <- function(
  mat_training,
  mat_test,
  fit,
  null_distributions,
  alpha = 0.01
) {

  k <- ncol(fit$W)

  # --- H scores ---
  H_train <- fit$H
  H_test  <- .predict_nmf(fit, mat_test)

  # --- p-values against null ---
  pval_train <- matrix(
    .null_pvalue(as.vector(H_train), null_distributions$H, null_distributions$probs),
    nrow = k
  )
  pval_test <- matrix(
    .null_pvalue(as.vector(H_test), null_distributions$H, null_distributions$probs),
    nrow = k
  )

  # --- BH correction per factor across cells ---
  padj_train <- apply(pval_train, 1, p.adjust, method = "BH")  # cells x factors
  padj_test  <- apply(pval_test,  1, p.adjust, method = "BH")

  # transpose back to k x cells
  padj_train <- t(padj_train)
  padj_test  <- t(padj_test)

  # --- binary membership ---
  membership_train <- pval_train < alpha
  membership_test  <- pval_test  < alpha


  # --- W gene membership ---
  n_genes <- nrow(mat_training) / 2

  pval_W_pos <- matrix(
    .null_pvalue(as.vector(fit$W[seq_len(n_genes), ]), null_distributions$W_pos, null_distributions$probs),
    nrow = n_genes
  )
  pval_W_neg <- matrix(
    .null_pvalue(as.vector(fit$W[seq(n_genes + 1, nrow(mat_training)), ]), null_distributions$W_neg, null_distributions$probs),
    nrow = n_genes
  )

  padj_W_pos <- apply(pval_W_pos, 2, p.adjust, method = "BH")  # genes x factors
  padj_W_neg <- apply(pval_W_neg, 2, p.adjust, method = "BH")

  gene_membership_pos <- pval_W_pos < alpha
  gene_membership_neg <- pval_W_neg < alpha
  
  out <- structure(
    list(
      membership_train = membership_train,
      membership_test  = membership_test,
      gene_pos         = gene_membership_pos,
      gene_neg         = gene_membership_neg,
      padj_train       = padj_train,
      padj_test        = padj_test,
      H_train          = H_train,
      H_test           = H_test,
      alpha            = alpha,
      padj_W_pos = padj_W_pos,
      padj_W_neg = padj_W_neg,
      pval_W_pos = pval_W_pos
    ),
    class = "bicluster_extracted"
  )
  out
}



.tail_enrichment <- function(observed, null_tails_list, threshold_prob, tail_threshold) {

  if (threshold_prob <= tail_threshold)
    stop(sprintf("threshold_prob %.4f must be above tail_threshold %.4f",
                 threshold_prob, tail_threshold))

  rescaled_prob <- (threshold_prob - tail_threshold) / (1 - tail_threshold)
  cutoff        <- quantile(unlist(null_tails_list), rescaled_prob)

  # number of exceedances per permutation
  null_exc_per_perm <- sapply(null_tails_list, function(tail) sum(tail > cutoff))
  
  null_exc_mean <- mean(null_exc_per_perm)
  null_exc_var  <- var(null_exc_per_perm)
  
  # observed number of exceedances
  obs_exc <- sum(observed > cutoff)
  
  list(
    cutoff       = cutoff,
    obs_exc      = obs_exc,
    null_exc     = null_exc_mean,
    null_exc_var = null_exc_var,
    enrichment   = obs_exc / null_exc_mean
  )
}

.find_cutoff_fdr <- function(observed, null_tails_list, tail_threshold, 
                               fdr_threshold = 0.2,
                               n_grid = 100) {
  
  # sweep over thresholds in the tail
  cutoffs <- quantile(unlist(null_tails_list), 
                      seq(0, 1, length.out = n_grid))

  results <- sapply(cutoffs, function(cut) {
    obs_exc       <- sum(observed > cut)
    null_exc_perm <- sapply(null_tails_list, function(tail) sum(tail > cut))
    null_exc_mean <- mean(null_exc_perm)
    null_exc_var  <- var(null_exc_perm)
    
    # noise fraction = expected null exceedances / observed exceedances
    if (obs_exc == 0) return(1)
    (null_exc_mean ) / obs_exc 
  })
  # we recruit a 95% confidence bound if we  do  (null_exc_mean + 1.96*sqrt(null_exc_var)) / obs_exc 
  # where we are sure that we are below fdr_threshold at 95% confidence!

  # find lowest cutoff where noise fraction < fdr_threshold
  below <- which(results < fdr_threshold)

  list(
    cutoff     = if (length(below) == 0) Inf else min(cutoffs[below]),
    noise_frac = data.frame(cutoff = cutoffs, noise_fraction = results)
  )
}


.null_exceedance <- function(null_density, threshold) {
  # integrate null density above threshold
  x <- null_density$x
  y <- null_density$y
  # only integrate above threshold
  above <- x >= threshold
  if (sum(above) < 2) return(0)
  # trapezoidal integration
  sum(diff(x[above]) * (y[above][-1] + y[above][-sum(above)]) / 2)
}


extract_biclusters_explicit <- function(
  mat_training,
  mat_test,
  fit,
  thrsh = 0.95
) {
  k       <- ncol(fit$W)
  n_genes <- nrow(fit$W) / 2

  # --- H scores ---
  H_train <- fit$H
  H_test  <- .predict_nmf(fit, mat_test)

  # --- W scores ---
  W_pos <- fit$W[seq_len(n_genes), ]
  W_neg <- fit$W[seq(n_genes + 1, 2 * n_genes), ]

  # --- cutoffs from observed quantile ---
  cutoff_W_pos <- quantile(as.vector(W_pos), thrsh)
  cutoff_W_neg <- quantile(as.vector(W_neg), thrsh)
  cutoff_H     <- quantile(as.vector(H_train), thrsh)

  # --- apply cutoff per factor ---
  gene_membership_pos <- W_pos > cutoff_W_pos
  gene_membership_neg <- W_neg > cutoff_W_neg
  membership_train    <- H_train > cutoff_H
  membership_test     <- H_test  > cutoff_H

  structure(
    list(
      membership_train  = membership_train,
      membership_test   = membership_test,
      gene_pos          = gene_membership_pos,
      gene_neg          = gene_membership_neg,
      H_train           = H_train,
      H_test            = H_test,
      cutoff_H          = cutoff_H,
      cutoff_W_pos      = cutoff_W_pos,
      cutoff_W_neg      = cutoff_W_neg,
      thrsh             = thrsh
    ),
    class = "bicluster_extracted_explicit"
  )
}


diagnose_biclusters <- function(bicluster_extracted) {
  if (!inherits(bicluster_extracted, "bicluster_extracted_explicit"))
    stop("'bicluster_extracted' must be the output of extract_biclusters()")

  # --- unpack ---
  k               <- nrow(bicluster_extracted$H_train)
  H_test          <- bicluster_extracted$H_test
  H_train         <- bicluster_extracted$H_train
  n_test          <- ncol(H_test)
  n_train         <- ncol(H_train)
  n_genes         <- nrow(bicluster_extracted$gene_pos)
  membership_test <- bicluster_extracted$membership_test
  fdr             <- bicluster_extracted$fdr


  # --- per factor summary ---
  summary_df <- data.frame(
    factor           = seq_len(k),
    n_cells_train    = rowSums(bicluster_extracted$membership_train),
    n_cells_test     = rowSums(bicluster_extracted$membership_test),
    frac_cells_train = rowSums(bicluster_extracted$membership_train)/n_train,
    frac_cells_test  = rowSums(bicluster_extracted$membership_test)/n_test,
    n_genes_pos      = colSums(bicluster_extracted$gene_pos),
    n_genes_neg      = colSums(bicluster_extracted$gene_neg),
    frac_genes_pos   = colSums(bicluster_extracted$gene_pos)/n_genes,
    frac_genes_neg   = colSums(bicluster_extracted$gene_neg)/n_genes
  )
  summary_df
}


.normalize_thr <- function(x, k) {
  if (length(x) == 1) {
    rep(x, k)
  } else if (length(x) == k) {
    x
  } else {
    stop(sprintf("must have length 1 or %d (got %d)", k, length(x)))
  }
}

extract_biclusters <- function(
  mat_training,
  mat_test,
  fit,
  thrsh_W_pos,
  thrsh_W_neg,
  thrsh_H
) {
  k       <- ncol(fit$W)
  n_genes <- nrow(fit$W) / 2

  thrsh_W_pos <- .normalize_thr(thrsh_W_pos, k)
  thrsh_W_neg <- .normalize_thr(thrsh_W_neg, k)
  thrsh_H     <- .normalize_thr(thrsh_H,     k)

  # --- H scores ---
  H_train <- fit$H
  H_test  <- .predict_nmf(fit, mat_test)

   # --- W scores ---
  W_pos <- fit$W[seq_len(n_genes), , drop = FALSE]
  W_neg <- fit$W[seq(n_genes + 1, 2 * n_genes), , drop = FALSE]

  cutoff_W_pos <- sapply(seq_len(ncol(W_pos)), function(j) {
    quantile(W_pos[, j], thrsh_W_pos[j], na.rm = TRUE)
  })

  cutoff_W_neg <- sapply(seq_len(ncol(W_neg)), function(j) {
    quantile(W_neg[, j], thrsh_W_neg[j], na.rm = TRUE)
  })

  cutoff_H <- sapply(seq_len(nrow(H_train)), function(j) {
    quantile(H_train[j, ], thrsh_H[j], na.rm = TRUE)
  })

  # --- apply cutoff per factor ---
  gene_membership_pos <- W_pos > matrix(cutoff_W_pos, nrow = nrow(W_pos), ncol = k, byrow = TRUE)
  gene_membership_neg <- W_neg > matrix(cutoff_W_neg, nrow = nrow(W_neg), ncol = k, byrow = TRUE)
  membership_train <- sweep(H_train, 1, cutoff_H, FUN = ">")
  membership_test  <- sweep(H_test,  1, cutoff_H, FUN = ">")

  rownames(gene_membership_pos) <- sub("^up_", "", rownames(gene_membership_pos))
  rownames(gene_membership_neg) <- sub("^down_", "", rownames(gene_membership_neg))

  list(
    gene_membership_pos = gene_membership_pos,
    gene_membership_neg = gene_membership_neg,
    membership_train = membership_train,
    membership_test = membership_test,
    cutoff_W_pos = cutoff_W_pos,
    cutoff_W_neg = cutoff_W_neg,
    cutoff_H = cutoff_H
  )
}


recruit_helper <- function(mat,
                           assay,
                           fit,
                           cutoffs = seq(0.99, 0.8, by = -0.01),
                           nnull = 20,
                           factor_id = 1,
                           seed = NULL,
                           use_irlba = TRUE,
                           work = 10,
                           maxit = 20,
                           tol = 1e-06) {
  if (!is.null(seed)) set.seed(seed)

  safe_z <- function(obs, mu, sd) {
    if (is.na(obs) || is.na(mu) || is.na(sd) || sd == 0) return(NA_real_)
    (obs - mu) / sd
  }

  emp_p_greater <- function(obs, null) {
    null <- null[is.finite(null)]
    if (!is.finite(obs) || length(null) == 0) return(NA_real_)
    (1 + sum(null >= obs, na.rm = TRUE)) / (length(null) + 1)
  }

  # Compute top 2 singular values and rank-1 fraction.
  # obs_pc1 here means sigma1^2 / ||X||_F^2, as before.
  rank2_score_fast <- function(X, v0 = NULL) {
    X <- as.matrix(X)
    if (anyNA(X)) X[is.na(X)] <- 0

    g <- nrow(X)
    c <- ncol(X)

    fro2 <- sum(X * X)
    if (!is.finite(fro2) || fro2 <= 0 || g < 2 || c < 2) {
      return(list(
        sigma1 = NA_real_,
        sigma2 = NA_real_,
        ratio1 = NA_real_,
        gap = NA_real_,
        u = NULL,
        v = NULL
      ))
    }

    # For tiny matrices, base svd is fine and often safer.
    if (min(g, c) <= 3) {
      s <- svd(X, nu = min(2, g), nv = min(2, c))
      d <- s$d
      sigma1 <- if (length(d) >= 1) abs(d[1]) else NA_real_
      sigma2 <- if (length(d) >= 2) abs(d[2]) else NA_real_
      ratio1 <- if (is.finite(sigma1)) (sigma1^2) / fro2 else NA_real_
      gap <- if (is.finite(sigma1) && is.finite(sigma2) && sigma2 > 0) sigma1 / sigma2 else NA_real_

      return(list(
        sigma1 = sigma1,
        sigma2 = sigma2,
        ratio1 = ratio1,
        gap = gap,
        u = if (!is.null(s$u)) s$u[, 1] else NULL,
        v = if (!is.null(s$v)) s$v[, 1] else NULL
      ))
    }

    # Preferred path: irlba if available.
    if (use_irlba && requireNamespace("irlba", quietly = TRUE)) {
      fit2 <- tryCatch({
        if (!is.null(v0) && length(v0) == c && all(is.finite(v0))) {
          irlba::irlba(X, nv = 2, nu = 2, v = v0, work = work, maxit = maxit, tol = tol)
        } else {
          irlba::irlba(X, nv = 2, nu = 2, work = work, maxit = maxit, tol = tol)
        }
      }, error = function(e) NULL)

      if (!is.null(fit2) && length(fit2$d) >= 2) {
        sigma1 <- abs(fit2$d[1])
        sigma2 <- abs(fit2$d[2])
        ratio1 <- (sigma1^2) / fro2
        gap <- if (is.finite(sigma2) && sigma2 > 0) sigma1 / sigma2 else NA_real_

        return(list(
          sigma1 = sigma1,
          sigma2 = sigma2,
          ratio1 = ratio1,
          gap = gap,
          u = fit2$u[, 1],
          v = fit2$v[, 1]
        ))
      }
    }

    # Fallback: simple power iteration for sigma1 only.
    # We still return sigma2 as NA in this fallback.
    v <- if (!is.null(v0) && length(v0) == c && all(is.finite(v0))) v0 else rnorm(c)
    v <- v / sqrt(sum(v^2))
    u <- rep(0, g)
    sigma_old <- NA_real_

    for (it in seq_len(10)) {
      u <- drop(X %*% v)
      unorm <- sqrt(sum(u^2))
      if (!is.finite(unorm) || unorm == 0) break
      u <- u / unorm

      v <- drop(crossprod(X, u))
      vnorm <- sqrt(sum(v^2))
      if (!is.finite(vnorm) || vnorm == 0) break
      v <- v / vnorm

      sigma <- abs(drop(crossprod(u, X %*% v)))
      if (is.finite(sigma_old) && abs(sigma - sigma_old) <= tol * max(1, abs(sigma_old))) {
        break
      }
      sigma_old <- sigma
    }

    sigma1 <- abs(drop(crossprod(u, X %*% v)))
    ratio1 <- (sigma1^2) / fro2

    if (sum(v, na.rm = TRUE) < 0) {
      u <- -u
      v <- -v
    }

    list(
      sigma1 = sigma1,
      sigma2 = NA_real_,
      ratio1 = ratio1,
      gap = NA_real_,
      u = u,
      v = v
    )
  }

  out <- vector("list", length(cutoffs))
  v_start <- NULL

  assay_rn <- rownames(assay)
  assay_cn <- colnames(assay)
  if (is.null(assay_rn) || is.null(assay_cn)) {
    stop("assay must have rownames and colnames.")
  }

  for (ii in seq_along(cutoffs)) {
    q <- cutoffs[ii]

    ext <- extract_biclusters_explicit(mat, mat, fit, q)
    rownames(ext$gene_pos) <- sub("^up_", "", rownames(ext$gene_pos))
    rownames(ext$gene_neg) <- sub("^down_", "", rownames(ext$gene_neg))

    genes <- union(
      names(which(ext$gene_pos[, factor_id] == TRUE)),
      names(which(ext$gene_neg[, factor_id] == TRUE))
    )
    cells <- names(which(ext$membership_train[factor_id, ] == TRUE))

    ng <- length(genes)
    ns <- length(cells)

    if (ng < 2 || ns < 2) {
      out[[ii]] <- data.frame(
        cutoff = q,
        ngenes = ng,
        nsamp = ns,
        obs_sigma1 = NA_real_,
        obs_sigma2 = NA_real_,
        obs_pc1 = NA_real_,
        obs_gap = NA_real_,
        null_pc1_mean = NA_real_,
        null_pc1_sd = NA_real_,
        null_gap_mean = NA_real_,
        null_gap_sd = NA_real_,
        p_pc1 = NA_real_,
        p_gap = NA_real_,
        z_pc1 = NA_real_,
        z_gap = NA_real_,
        row.names = NULL
      )
      next
    }

    gidx <- match(genes, assay_rn)
    cidx <- match(cells, assay_cn)
    gidx <- gidx[!is.na(gidx)]
    cidx <- cidx[!is.na(cidx)]

    if ( (length(genes_pos) + length(genes_neg)) < 2 || length(cidx) < 2) {
      out[[ii]] <- data.frame(
        cutoff = q,
        ngenes = length(gidx),
        nsamp = length(cidx),
        obs_sigma1 = NA_real_,
        obs_sigma2 = NA_real_,
        obs_pc1 = NA_real_,
        obs_gap = NA_real_,
        null_pc1_mean = NA_real_,
        null_pc1_sd = NA_real_,
        null_gap_mean = NA_real_,
        null_gap_sd = NA_real_,
        p_pc1 = NA_real_,
        p_gap = NA_real_,
        z_pc1 = NA_real_,
        z_gap = NA_real_,
        row.names = NULL
      )
      next
    }

    Xobs <- assay[gidx, cidx, drop = FALSE]
    #Xobs <- sweep(Xobs, 1, rowMeans(Xobs), "-") # Takes care of the sweep also for later!

    obs_fit <- rank2_score_fast(Xobs, v0 = if (!is.null(v_start) && length(v_start) == ncol(Xobs)) v_start else NULL)
    obs_sigma1 <- obs_fit$sigma1
    obs_sigma2 <- obs_fit$sigma2
    obs_pc1 <- obs_fit$ratio1
    obs_gap <- obs_fit$gap

    v_start <- obs_fit$v

    null_pc1 <- numeric(nnull)
    null_gap <- numeric(nnull)

    for (jj in seq_len(nnull)) {
      gnull <- sample.int(nrow(assay), length(gidx))
      snull <- sample.int(ncol(assay), length(cidx))
      Xnull <- assay[gnull, snull, drop = FALSE]
      
      #Xnull <- Xobs
      #for (rr in seq_len(nrow(Xnull))) {
      #  Xnull[rr, ] <- sample(Xnull[rr, ])
      #}

      nf <- rank2_score_fast(Xnull, v0 = NULL)
      null_pc1[jj] <- nf$ratio1
      null_gap[jj] <- nf$gap
    }

    m_pc1 <- mean(null_pc1, na.rm = TRUE)
    s_pc1 <- sd(null_pc1, na.rm = TRUE)
    m_gap <- mean(null_gap, na.rm = TRUE)
    s_gap <- sd(null_gap, na.rm = TRUE)

    out[[ii]] <- data.frame(
      cutoff = q,
      ngenes = length(gidx),
      nsamp = length(cidx),

      obs_sigma1 = obs_sigma1,
      obs_sigma2 = obs_sigma2,
      obs_pc1 = obs_pc1,
      obs_gap = obs_gap,

      null_pc1_mean = m_pc1,
      null_pc1_sd = s_pc1,
      null_gap_mean = m_gap,
      null_gap_sd = s_gap,

      p_pc1 = emp_p_greater(obs_pc1, null_pc1),
      p_gap = emp_p_greater(obs_gap, null_gap),

      z_pc1 = safe_z(obs_pc1, m_pc1, s_pc1),
      z_gap = safe_z(obs_gap, m_gap, s_gap),

      row.names = NULL
    )
  }

  do.call(rbind, out)
}





recruit_helper_simpler <- function(mat,
                           assay,
                           fit,
                           cutoffs = seq(0.99, 0.8, by = -0.01),
                           nnull = 20,
                           factor_id = 1,
                           seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  score_summary <- function(X) {
    gene_vec <- abs(rowMeans(X, na.rm = TRUE))
    cell_vec <- abs(colMeans(abs(X), na.rm = TRUE))

    data.frame(
      mean_gene = mean(gene_vec, na.rm = TRUE),
      mean_cell = mean(cell_vec, na.rm = TRUE),
      var_gene  = if (length(gene_vec) > 1) var(gene_vec, na.rm = TRUE) else NA_real_,
      var_cell  = if (length(cell_vec) > 1) var(cell_vec, na.rm = TRUE) else NA_real_
    )
  }

  safe_z <- function(obs, mu, sd) {
    if (is.na(obs) || is.na(mu) || is.na(sd) || sd == 0) return(NA_real_)
    (obs - mu) / sd
  }

  out <- vector("list", length(cutoffs))

  for (ii in seq_along(cutoffs)) {
    q <- cutoffs[ii]

    ext <- extract_biclusters_explicit(mat, mat, fit, q)
    rownames(ext$gene_pos) <- sub("^up_", "", rownames(ext$gene_pos))
    rownames(ext$gene_neg) <- sub("^down_", "", rownames(ext$gene_neg))

    genes <- union(
      names(which(ext$gene_pos[, factor_id] == TRUE)),
      names(which(ext$gene_neg[, factor_id] == TRUE))
    )
    cells <- names(which(ext$membership_train[factor_id, ] == TRUE))

    ng <- length(genes)
    ns <- length(cells)

    if (ng < 2 || ns < 2) {
      out[[ii]] <- data.frame(
        cutoff = q,
        ngenes = ng,
        nsamp = ns,
        obs_mean_gene = NA_real_,
        obs_mean_cell = NA_real_,
        obs_var_gene = NA_real_,
        obs_var_cell = NA_real_,
        null_mean_gene_mean = NA_real_,
        null_mean_gene_sd = NA_real_,
        null_mean_cell_mean = NA_real_,
        null_mean_cell_sd = NA_real_,
        null_var_gene_mean = NA_real_,
        null_var_gene_sd = NA_real_,
        null_var_cell_mean = NA_real_,
        null_var_cell_sd = NA_real_,
        z_mean_gene = NA_real_,
        z_mean_cell = NA_real_,
        z_var_gene = NA_real_,
        z_var_cell = NA_real_,
        row.names = NULL
      )
      next
    }

    Xobs <- assay[genes, cells, drop = FALSE]
    obs <- score_summary(Xobs)

    null_stats <- replicate(nnull, {
      gnull <- sample(seq_len(nrow(assay)), ng)
      snull <- sample(seq_len(ncol(assay)), ns)
      Xnull <- assay[gnull, snull, drop = FALSE]
      score_summary(Xnull)
    }, simplify = FALSE)

    null_df <- do.call(rbind, null_stats)

    m_mean_gene <- mean(null_df$mean_gene, na.rm = TRUE)
    s_mean_gene <- sd(null_df$mean_gene, na.rm = TRUE)

    m_mean_cell <- mean(null_df$mean_cell, na.rm = TRUE)
    s_mean_cell <- sd(null_df$mean_cell, na.rm = TRUE)

    m_var_gene <- mean(null_df$var_gene, na.rm = TRUE)
    s_var_gene <- sd(null_df$var_gene, na.rm = TRUE)

    m_var_cell <- mean(null_df$var_cell, na.rm = TRUE)
    s_var_cell <- sd(null_df$var_cell, na.rm = TRUE)

    out[[ii]] <- data.frame(
      cutoff = q,
      ngenes = ng,
      nsamp = ns,

      obs_mean_gene = obs$mean_gene,
      obs_mean_cell = obs$mean_cell,
      obs_var_gene  = obs$var_gene,
      obs_var_cell  = obs$var_cell,

      null_mean_gene_mean = m_mean_gene,
      null_mean_gene_sd = s_mean_gene,
      null_mean_cell_mean = m_mean_cell,
      null_mean_cell_sd = s_mean_cell,
      null_var_gene_mean = m_var_gene,
      null_var_gene_sd = s_var_gene,
      null_var_cell_mean = m_var_cell,
      null_var_cell_sd = s_var_cell,

      z_mean_gene = safe_z(obs$mean_gene, m_mean_gene, s_mean_gene),
      z_mean_cell = safe_z(obs$mean_cell, m_mean_cell, s_mean_cell),
      z_var_gene  = safe_z(obs$var_gene,  m_var_gene,  s_var_gene),
      z_var_cell  = safe_z(obs$var_cell,  m_var_cell,  s_var_cell),

      row.names = NULL
    )
  }

  do.call(rbind, out)
}




.crit <- function(df) {
  with(df, z_mean > 3 & z_pc1 > 3 & (z_mean + z_pc1) > 10)
}

loop_factors <- function(mat,
                         assay,
                         fit,
                         cutoffs = seq(0.99, 0.8, by = -0.01),
                         nnull = 10,
                         factor_ids = NULL,
                         seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Guess number of factors if not provided
  if (is.null(factor_ids)) {
    factor_ids <- seq_len(nrow(fit$H))
  }

  res <- vector("list", length(factor_ids))

  for (jj in seq_along(factor_ids)) {
    f <- factor_ids[jj]

    stats <- recruit_helper(
      mat = mat,
      assay = assay,
      fit = fit,
      cutoffs = cutoffs,
      nnull = nnull,
      factor_id = f,
      seed = seed
    )

    ok <- .crit(stats)

    if (length(ok) != nrow(stats)) {
      stop(" must return a logical vector with one value per cutoff row.")
    }

    if (any(ok, na.rm = TRUE)) {
      # last cutoff in the cutoff sequence that still fulfills the criterion
      idx <- tail(which(ok), 1)
      chosen <- stats[idx, , drop = FALSE]
      chosen$factor_id <- f
    } else {
      chosen <- data.frame(
        factor_id = f,
        cutoff = NA_real_,
        ngenes = NA_integer_,
        nsamp = NA_integer_,
        obs_mean = NA_real_,
        obs_pc1 = NA_real_,
        null_mean_mean = NA_real_,
        null_mean_sd = NA_real_,
        null_pc1_mean = NA_real_,
        null_pc1_sd = NA_real_,
        z_mean = NA_real_,
        z_pc1 = NA_real_
      )
    }

    res[[jj]] <- chosen
  }

  do.call(rbind, res)
}



loop_factors_overview <- function(mat,
                                  assay,
                                  fit,
                                  cutoffs = seq(0.99, 0.8, by = -0.01),
                                  nnull = 10,
                                  factor_ids = NULL,
                                  seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  if (is.null(factor_ids)) {
    factor_ids <- seq_len(nrow(fit$H))
  }

  res <- vector("list", length(factor_ids))

  for (jj in seq_along(factor_ids)) {
    f <- factor_ids[jj]

    stats <- recruit_helper_simpler(
      mat = mat,
      assay = assay,
      fit = fit,
      cutoffs = cutoffs,
      nnull = nnull,
      factor_id = f,
      seed = seed
    )

    stats$factor_id <- f
    stats$cutoff_rank <- seq_len(nrow(stats))

    res[[jj]] <- stats
  }

  out <- do.call(rbind, res)
  rownames(out) <- NULL

  out <- out[order(out$factor_id, -out$cutoff), ]
  rownames(out) <- NULL

  out
}