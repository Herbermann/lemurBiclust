# ============================================================
# HELPERS
# ============================================================

.emp_p_greater <- function(obs, null) {
  null <- null[is.finite(null)]

  if (!is.finite(obs) || length(null) == 0) {
    return(NA_real_)
  }

  (1 + sum(null >= obs, na.rm = TRUE)) /
    (length(null) + 1)
}

.safe_z <- function(obs, mu, sd) {
  if (is.na(obs) || is.na(mu) || is.na(sd) || sd == 0) {
    return(NA_real_)
  }

  (obs - mu) / sd
}

.fast_extract <- function(ranked_list, q) {
  n <- length(ranked_list)
  k <- ceiling((1 - q) * n)
  if (k < 1) {
    return(character(0))
  }
  names(ranked_list)[seq_len(k)]
}

.scalar_or_na <- function(x) {
  if (is.null(x) || length(x) == 0) {
    NA_real_
  } else {
    x[1]
  }
}

.score_bicluster <- function(
    X,
    v0 = NULL,
    use_irlba = TRUE,
    return_gap = FALSE,
    work = 10,
    maxit = 20,
    tol = 1e-06
) {

  X <- as.matrix(X)

  if (anyNA(X)) {
    X[is.na(X)] <- 0
  }

  g <- nrow(X)
  c <- ncol(X)

  fro2 <- sum(X * X)

  if (!is.finite(fro2) ||
      fro2 <= 0 ||
      g < 2 ||
      c < 2) {

    return(list(
      sigma1 = NA_real_,
      sigma2 = NA_real_,
      ratio1 = NA_real_,
      gap = NA_real_,
      v = NULL
    ))
  }

  n_eigenvectors <- if (return_gap) 2 else 1

  if (use_irlba &&
      requireNamespace("irlba", quietly = TRUE)) {

    fit2 <- tryCatch({

      if (!is.null(v0) &&
          length(v0) == c &&
          all(is.finite(v0))) {

        irlba::irlba(
          X,
          nv = n_eigenvectors,
          nu = 0,
          v = v0,
          work = work,
          maxit = maxit,
          tol = tol
        )

      } else {

        irlba::irlba(
          X,
          nv = n_eigenvectors,
          nu = 0,
          work = work,
          maxit = maxit,
          tol = tol
        )
      }

      }, error = function(e) {
        print(e)
        NULL
      })

    if (!is.null(fit2) &&
        length(fit2$d) >= 1) {

      sigma1 <- abs(fit2$d[1])

      sigma2 <- if (return_gap &&
                    length(fit2$d) >= 1) {
        abs(fit2$d[2])
      } else {
        NA_real_
      }

      gap <- if (return_gap &&
                 is.finite(sigma2) &&
                 sigma2 > 0) {
        sigma1 / sigma2
      } else {
        NA_real_
      }

      v_out <- if (is.matrix(fit2$v)) {
        fit2$v[, 1]
      } else {
        fit2$v
      }

      return(list(
        sigma1 = sigma1,
        sigma2 = sigma2,
        ratio1 = (sigma1^2) / fro2,
        gap = gap,
        v = v_out
      ))
    }
  }

  return(list(
    sigma1 = NA_real_,
    sigma2 = NA_real_,
    ratio1 = NA_real_,
    gap = NA_real_,
    v = NULL
  ))
}



.null_manager <- function(
    assay,
    nnull = 100,
    use_irlba = TRUE,
    work = 10,
    maxit = 20,
    tol = 1e-06,
    bin_genes = 25,
    bin_cells = 100
) {

  # ------------------------------------------------------------
  # INTERNAL CACHE
  # ------------------------------------------------------------

  cache <- new.env(parent = emptyenv())

  # ------------------------------------------------------------
  # DIMENSION BINNING
  # ------------------------------------------------------------

  make_key <- function(ng, nc) {

    ng_bin <- max(
      bin_genes,
      round(ng / bin_genes) * bin_genes
    )   
      
    nc_bin <- max(
      bin_cells,
      round(nc / bin_cells) * bin_cells
    )
    
    paste0(ng_bin, "x", nc_bin)
  }

  # ------------------------------------------------------------
  # BUILD NULL DISTRIBUTION
  # ------------------------------------------------------------

  compute_null <- function(ng, nc, return_gap = FALSE) {

    null_pc1 <- numeric(nnull)

    if (return_gap) {
      null_gap <- numeric(nnull)
    } else {
      null_gap <- NULL
    }

    for (jj in seq_len(nnull)) {

      gnull <- sample.int(nrow(assay), ng)
      snull <- sample.int(ncol(assay), nc)

      Xnull <- assay[gnull, snull, drop = FALSE]

      nf <- .score_bicluster(
        Xnull,
        use_irlba = use_irlba,
        return_gap = return_gap,
        work = work,
        maxit = maxit,
        tol = tol
      )

      null_pc1[jj] <- nf$ratio1

      if (return_gap) {
        null_gap[jj] <- nf$gap
      }
    }

    out <- list(

      pc1 = null_pc1,

      pc1_mean = mean(null_pc1, na.rm = TRUE),
      pc1_sd = sd(null_pc1, na.rm = TRUE)

    )

    if (return_gap) {

      out$gap <- null_gap

      out$gap_mean <- mean(null_gap, na.rm = TRUE)
      out$gap_sd <- sd(null_gap, na.rm = TRUE)
    }

    out
  }

  # ------------------------------------------------------------
  # PUBLIC GET METHOD
  # ------------------------------------------------------------

  get_null <- function(ng, nc, return_gap = FALSE) {

    key <- paste0(
      make_key(ng, nc),
      "_gap_",
      return_gap
    )

    # CACHE HIT
    if (exists(key, envir = cache, inherits = FALSE)) {

      return(get(key, envir = cache))
    }

    # CACHE MISS
    out <- compute_null(
      ng = ng,
      nc = nc,
      return_gap = return_gap
    )

    assign(key, out, envir = cache)

    out
  }

  # ------------------------------------------------------------
  # OPTIONAL CACHE INSPECTION
  # ------------------------------------------------------------

  cache_keys <- function() {
    ls(envir = cache)
  }

  clear_cache <- function() {
    rm(list = ls(envir = cache), envir = cache)
  }

  # ------------------------------------------------------------
  # RETURN MANAGER OBJECT
  # ------------------------------------------------------------

  list(
    get_null = get_null,
    keys = cache_keys,
    clear = clear_cache
  )
}


.factor_recruiter <- function(
    mat,
    assay,
    ranked_W_pos,
    ranked_W_neg,
    ranked_H,
    cutoff_pos,
    cutoff_neg,
    cutoff_cell,
    factor_id,
    manager,
    return_gap = FALSE,
    v0 = NULL
){

  assay_rn <- rownames(assay)
  assay_cn <- colnames(assay)

  if (is.null(assay_rn) ||
      is.null(assay_cn)) {

    stop("assay needs rownames and colnames")
  }

  # Lookup table for names <-> indices
  gene_map <- setNames(seq_along(assay_rn), assay_rn)
  cell_map <- setNames(seq_along(assay_cn), assay_cn)
  
  # recruit genes and cells for this factor at this cutoff
  genes_pos <- .fast_extract(ranked_W_pos, cutoff_pos)
  genes_neg <- .fast_extract(ranked_W_neg, cutoff_neg)
  cells     <- .fast_extract(ranked_H, cutoff_cell)

  genes <- union(genes_pos, genes_neg)

  gidx <- unname(gene_map[genes])
  cidx <- unname(cell_map[cells])

  gidx <- gidx[!is.na(gidx)]
  cidx <- cidx[!is.na(cidx)]

  
  # Build at rate bicluster!
  Xobs <- assay[gidx, cidx, drop = FALSE]

  current <- .score_bicluster(
    Xobs,
    v0 = v0,
    return_gap
  )

  # Get results from null_manager for corresponding null!

  null_results <- manager$get_null(length(gidx), length(cidx))


   list(
      q_pos = cutoff_pos,
      q_neg = cutoff_neg,
      q_cell = cutoff_cell,

      ngenes_pos = length(genes_pos),
      ngenes_neg = length(genes_neg),

      ngenes_total = length(gidx),
      ncells = length(cidx),

      sigma1 = current$sigma1,
      sigma2 = current$sigma2,

      pc1 = current$ratio1,
      gap = current$gap,

      null_pc1_mean = null_results$pc1_mean,
      null_pc1_sd = null_results$pc1_sd,

      null_gap_mean = null_results$gap_mean,
      null_gap_sd = null_results$gap_sd,

      p_pc1 = .emp_p_greater(
        current$ratio1,
        null_results$pc1
      ),

      p_gap = .emp_p_greater(
        current$gap,
        null_results$gap
      ),
      
      v = current$v

  )
}




.global_recruiter <- function(
    mat,
    assay,
    fit,
    factors = NULL,

    cutoffs_pos = seq(0.99, 0.8, by = -0.01),
    cutoffs_neg = seq(0.99, 0.8, by = -0.01),
    cutoffs_cell = seq(0.99, 0.8, by = -0.01),

    cut_pc_abs = 0.1,
    alpha = 0.05,

    nnull = 100,
    use_irlba = TRUE,
    return_gap = FALSE
) {

  # ============================================================
  # NULL MANAGER
  # ============================================================

  manager <- .null_manager(
    assay = assay,
    nnull = nnull,
    use_irlba = use_irlba
  )

  # ============================================================
  # FACTOR SET
  # ============================================================

  if (is.null(factors)) {
    factors <- seq_len(nrow(fit$H))
  }

  # ============================================================
  # PRECOMPUTE RANKINGS
  # ============================================================

  rankings <- lapply(
    factors,
    function(k) {

      n_genes <- nrow(fit$W) / 2

      W_pos <- sort(
        fit$W[1:n_genes, k],
        decreasing = TRUE
      )

      W_neg <- sort(
        fit$W[(n_genes + 1):(2 * n_genes), k],
        decreasing = TRUE
      )

      H <- sort(
        fit$H[k, ],
        decreasing = TRUE
      )

      names(W_pos) <- sub("^up_", "", names(W_pos))
      names(W_neg) <- sub("^down_", "", names(W_neg))

      list(
        W_pos = W_pos,
        W_neg = W_neg,
        H = H
      )
    }
  )

  # ============================================================
  # RESULTS
  # ============================================================

  results <- vector("list", length(factors))

  # ============================================================
  # LOOP OVER FACTORS
  # ============================================================

  for (factor_idx in seq_along(factors)) {

    cur_factor <- factors[factor_idx]

    # ----------------------------------------------------------
    # INITIALIZE RECRUITMENT STATE
    # ----------------------------------------------------------

    idx_pos <- 1
    idx_neg <- 1
    idx_cell <- 1

    warm_v <- NULL

    accepted <- list()
    proposals <- list()

    idx_accepted <- 1
    idx_proposal <- 1

    step <- 1

    # ==========================================================
    # FACTOR TRAJECTORY LOOP
    # ==========================================================

    repeat {

      # --------------------------------------------------------
      # CURRENT FROZEN STATE
      # --------------------------------------------------------

      cur_pos <- idx_pos
      cur_neg <- idx_neg
      cur_cell <- idx_cell

      moved <- FALSE

      # --------------------------------------------------------
      # EVALUATE CURRENT ACCEPTED STATE
      # --------------------------------------------------------

      current <- .factor_recruiter(
        mat = mat,
        assay = assay,

        ranked_W_pos = rankings[[factor_idx]]$W_pos,
        ranked_W_neg = rankings[[factor_idx]]$W_neg,
        ranked_H = rankings[[factor_idx]]$H,

        cutoff_pos = cutoffs_pos[cur_pos],
        cutoff_neg = cutoffs_neg[cur_neg],
        cutoff_cell = cutoffs_cell[cur_cell],

        factor_id = cur_factor,

        manager = manager,

        return_gap = return_gap,
        v0 = warm_v
      )

      # failed recruitment state
      if (is.null(current)) {
        break
      }

      # update warm start
      warm_v <- current$v


      # get the null
      ng <- length(current$ngenes_total)
      nc <- length(current$ncells)
      null_results <- manager$get_null(ng, nc, return_gap = TRUE)


      # --------------------------------------------------------
      # STORE ACCEPTED CURRENT STATE
      # --------------------------------------------------------

      # Format for safety
      sigma2 <- .scalar_or_na(current$sigma2)
      gap    <- .scalar_or_na(current$gap)

      p_gap  <- .scalar_or_na(
        .emp_p_greater(
          current$gap,
          null_results$gap
        )
      )

      z_gap <- .scalar_or_na(
        .safe_z(
          current$gap,
          null_results$gap_mean,
          null_results$gap_sd
        )
      )
      
      # Safety net.
      sigma1 <- .scalar_or_na(current$sigma1)
      sigma2 <- .scalar_or_na(current$sigma2)

      pc1 <- .scalar_or_na(current$pc1)
      gap <- .scalar_or_na(current$gap)

      p_pc1 <- .scalar_or_na(current$p_pc1)
      p_gap <- .scalar_or_na(current$p_gap)

      z_pc1 <- .scalar_or_na(current$z_pc1)
      z_gap <- .scalar_or_na(current$z_gap)


      accepted[[idx_accepted]] <- data.frame(

        factor = cur_factor,
        step = step,

        cutoff_pos = current$q_pos,
        cutoff_neg = current$q_neg,
        cutoff_cell = current$q_cell,

        ngenes_pos = current$ngenes_pos,
        ngenes_neg = current$ngenes_neg,

        ngenes_total = current$ngenes_total,
        ncells = current$ncells,

        sigma1 = sigma1,
        sigma2 = sigma2,

        pc1 = pc1,
        gap = gap,

        p_pc1 = p_pc1,
        p_gap = p_gap,

        z_pc1 = z_pc1,
        z_gap = z_gap,

        row.names = NULL
      )

      idx_accepted <- idx_accepted + 1

      # ========================================================
      # PROPOSAL 1: POSITIVE GENES
      # ========================================================

      if (cur_pos < length(cutoffs_pos)) {

        proposal <- .factor_recruiter(
          mat = mat,
          assay = assay,

          ranked_W_pos = rankings[[factor_idx]]$W_pos,
          ranked_W_neg = rankings[[factor_idx]]$W_neg,
          ranked_H = rankings[[factor_idx]]$H,

          cutoff_pos = cutoffs_pos[cur_pos + 1],
          cutoff_neg = cutoffs_neg[cur_neg],
          cutoff_cell = cutoffs_cell[cur_cell],

          factor_id = cur_factor,

          manager = manager,

          return_gap = return_gap,
          v0 = warm_v
        )

        proposal_accept <-
          !is.null(proposal) &&
          proposal$p_pc1 < alpha &&
          proposal$pc1 > cut_pc_abs

        proposals[[idx_proposal]] <- data.frame(

          factor = cur_factor,
          step = step,

          axis = "pos",

          current_pos = cutoffs_pos[cur_pos],
          proposed_pos = cutoffs_pos[cur_pos + 1],

          current_neg = cutoffs_neg[cur_neg],
          proposed_neg = NA_real_,

          current_cell = cutoffs_cell[cur_cell],
          proposed_cell = NA_real_,

          accepted = proposal_accept,

          proposal_pc1 <- if (is.null(proposal)) {
              NA_real_
          } else {
              proposal$pc1
          },

          proposal_gap <- if (is.null(proposal)) {
              NA_real_
          } else {
              proposal$gap
          },

          proposal_p <- if (is.null(proposal)) {
              NA_real_
          } else {
              proposal$p_pc1
          },

          row.names = NULL
        )

        idx_proposal <- idx_proposal + 1

        if (proposal_accept) {
          idx_pos <- idx_pos + 1
          moved <- TRUE
        }
      }

      # ========================================================
      # PROPOSAL 2: NEGATIVE GENES
      # ========================================================

      if (cur_neg < length(cutoffs_neg)) {

        proposal <- .factor_recruiter(
          mat = mat,
          assay = assay,

          ranked_W_pos = rankings[[factor_idx]]$W_pos,
          ranked_W_neg = rankings[[factor_idx]]$W_neg,
          ranked_H = rankings[[factor_idx]]$H,

          cutoff_pos = cutoffs_pos[cur_pos],
          cutoff_neg = cutoffs_neg[cur_neg + 1],
          cutoff_cell = cutoffs_cell[cur_cell],

          factor_id = cur_factor,

          manager = manager,

          return_gap = return_gap,
          v0 = warm_v
        )

        proposal_accept <-
          !is.null(proposal) &&
          proposal$p_pc1 < alpha &&
          proposal$pc1 > cut_pc_abs

        proposals[[idx_proposal]] <- data.frame(

          factor = cur_factor,
          step = step,

          axis = "neg",

          current_pos = cutoffs_pos[cur_pos],
          proposed_pos = NA_real_,

          current_neg = cutoffs_neg[cur_neg],
          proposed_neg = cutoffs_neg[cur_neg + 1],

          current_cell = cutoffs_cell[cur_cell],
          proposed_cell = NA_real_,

          accepted = proposal_accept,

          proposal_pc1 <- if (is.null(proposal)) {
              NA_real_
          } else {
              proposal$pc1
          },

          proposal_gap <- if (is.null(proposal)) {
              NA_real_
          } else {
              proposal$gap
          },

          proposal_p <- if (is.null(proposal)) {
              NA_real_
          } else {
              proposal$p_pc1
          },

          row.names = NULL
        )

        idx_proposal <- idx_proposal + 1

        if (proposal_accept) {
          idx_neg <- idx_neg + 1
          moved <- TRUE
        }
      }

      # ========================================================
      # PROPOSAL 3: CELLS
      # ========================================================

      if (cur_cell < length(cutoffs_cell)) {

        proposal <- .factor_recruiter(
          mat = mat,
          assay = assay,

          ranked_W_pos = rankings[[factor_idx]]$W_pos,
          ranked_W_neg = rankings[[factor_idx]]$W_neg,
          ranked_H = rankings[[factor_idx]]$H,

          cutoff_pos = cutoffs_pos[cur_pos],
          cutoff_neg = cutoffs_neg[cur_neg],
          cutoff_cell = cutoffs_cell[cur_cell + 1],

          factor_id = cur_factor,

          manager = manager,

          return_gap = return_gap,
          v0 = warm_v
        )

        proposal_accept <-
          !is.null(proposal) &&
          proposal$p_pc1 < alpha &&
          proposal$pc1 > cut_pc_abs

        proposals[[idx_proposal]] <- data.frame(

          factor = cur_factor,
          step = step,

          axis = "cell",

          current_pos = cutoffs_pos[cur_pos],
          proposed_pos = NA_real_,

          current_neg = cutoffs_neg[cur_neg],
          proposed_neg = NA_real_,

          current_cell = cutoffs_cell[cur_cell],
          proposed_cell = cutoffs_cell[cur_cell + 1],

          accepted = proposal_accept,

          proposal_pc1 <- if (is.null(proposal)) {
              NA_real_
          } else {
              proposal$pc1
          },

          proposal_gap <- if (is.null(proposal)) {
              NA_real_
          } else {
              proposal$gap
          },

          proposal_p <- if (is.null(proposal)) {
              NA_real_
          } else {
              proposal$p_pc1
          },

          row.names = NULL
        )

        idx_proposal <- idx_proposal + 1

        if (proposal_accept) {
          idx_cell <- idx_cell + 1
          moved <- TRUE
        }
      }

      # --------------------------------------------------------
      # TERMINATION
      # --------------------------------------------------------

      if (!moved) {
        break
      }

      step <- step + 1
    }

    # ==========================================================
    # STORE FACTOR RESULT
    # ==========================================================

    results[[factor_idx]] <- list(

      factor = cur_factor,

      accepted = do.call(rbind, accepted),

      proposals = do.call(rbind, proposals)
    )
  }

  biclusters <- do.call(rbind, lapply(seq_along(results), function(i) {

    acc <- results[[i]]$accepted

    if (is.null(acc) || nrow(acc) == 0) {
      return(data.frame(
        factor = i,
        cutoff_pos = NA_real_,
        cutoff_neg = NA_real_,
        cutoff_cell = NA_real_,
        ngenes_pos = NA_real_,
        ngenes_neg = NA_real_,
        ncells = NA_real_,
        pc1 = NA_real_
      ))
    }

    cbind(
      factor = i,
      acc[
        nrow(acc),
        c(
          "cutoff_pos",
          "cutoff_neg",
          "cutoff_cell",
          "ngenes_pos",
          "ngenes_neg",
          "ncells",
          "pc1"
        ),
        drop = FALSE
      ]
    )
  }))

  out = list(
    recruitment_log = results,
    biclusters = biclusters
  )
  
}