.TRANSFORMS <- list(
  none  = function(...) identity,
  tanh  = function(w = 1, ...) function(x) tanh(w * x) / w,
  asinh = function(w = 1, ...) function(x) asinh(w * x) / w
)

.apply_fixed_noise_clip <- function(mat, thrsh) {
  if (thrsh <= 0) return(mat)

  if (inherits(mat, "sparseMatrix")) {
    out <- mat
    if (length(out@x)) {
      out@x[abs(out@x) < thrsh] <- 0
      out <- Matrix::drop0(out)
    }
    out
  } else {
    out <- mat
    out[abs(out) < thrsh] <- 0
    out
  }
}

.split_pos <- function(mat) {
  if (inherits(mat, "sparseMatrix")) {
    out <- mat
    if (length(out@x)) {
      out@x[out@x < 0] <- 0
      out <- Matrix::drop0(out)
    }
    out
  } else {
    pmax(mat, 0)
  }
}

.split_neg <- function(mat) {
  if (inherits(mat, "sparseMatrix")) {
    out <- mat
    if (length(out@x)) {
      out@x <- pmax(-out@x, 0)
      out <- Matrix::drop0(out)
    }
    out
  } else {
    pmax(-mat, 0)
  }
}

.abs_quantile <- function(mat, p) {
  if (inherits(mat, "sparseMatrix")) {
    n <- nrow(mat) * ncol(mat)
    nz <- abs(mat@x)
    nnz <- length(nz)
    n0 <- n - nnz

    if (n0 > 0 && p <= n0 / n) return(0)
    if (nnz == 0) return(0)

    p2 <- (p * n - n0) / nnz
    p2 <- min(max(p2, 0), 1)
    as.numeric(stats::quantile(nz, probs = p2, names = FALSE, na.rm = TRUE))
  } else {
    as.numeric(stats::quantile(abs(mat), probs = p, names = FALSE, na.rm = TRUE))
  }
}

.pos_quantile <- function(mat, p) {
  if (inherits(mat, "sparseMatrix")) {
    nz <- mat@x[mat@x > 0]
    if (!length(nz)) return(0)
    as.numeric(stats::quantile(nz, probs = p, names = FALSE, na.rm = TRUE))
  } else {
    nz <- mat[mat > 0]
    if (!length(nz)) return(0)
    as.numeric(stats::quantile(nz, probs = p, names = FALSE, na.rm = TRUE))
  }
}

.compute_weight <- function(mat_half, clip_max, clip_max_value) {
  if (clip_max == "none") return(1)

  if (is.null(clip_max_value) || length(clip_max_value) != 1 || !is.finite(clip_max_value)) {
    stop("'clip_max_value' must be a finite scalar when clip_max != 'none'")
  }

  if (clip_max == "fixed") {
    if (clip_max_value <= 0) stop("'clip_max_value' must be > 0 when clip_max = 'fixed'")
    return(1 / clip_max_value)
  }

  if (clip_max == "quantile") {
    q <- .pos_quantile(mat_half, clip_max_value)
    if (q <= 0) return(1)
    return(1 / q)
  }

  1
}

.apply_transform_sparse <- function(mat, transform, w) {
  if (transform == "none") return(mat)

  out <- mat
  if (length(out@x)) {
    if (transform == "tanh") {
      out@x <- tanh(w * out@x) / w
    } else if (transform == "asinh") {
      out@x <- asinh(w * out@x) / w
    } else {
      stop("unknown transform")
    }
  }
  out
}

.apply_transform_dense <- function(mat, transform, w) {
  if (transform == "none") return(mat)
  if (transform == "tanh")  return(tanh(w * mat) / w)
  if (transform == "asinh")  return(asinh(w * mat) / w)
  stop("unknown transform")
}


preprocess_assay_fit <- function(
  lemur_fit_class,
  use_assay,
  transform        = c("none", "tanh", "asinh"),
  clip_noise       = c("none", "fixed", "quantile"),
  clip_noise_value = NULL,
  clip_max         = c("none", "quantile", "fixed"),
  clip_max_value   = NULL
) {
  transform  <- match.arg(transform)
  clip_noise <- match.arg(clip_noise)
  clip_max   <- match.arg(clip_max)

  if (!use_assay %in% assayNames(lemur_fit_class)) {
    stop(sprintf(
      "'use_assay' must be one of: %s",
      paste(assayNames(lemur_fit_class), collapse = ", ")
    ))
  }

  mat_raw <- assay(lemur_fit_class, use_assay)

  noise_thr <- NULL
  if (clip_noise == "fixed") {
    if (is.null(clip_noise_value) || length(clip_noise_value) != 1 || !is.finite(clip_noise_value)) {
      stop("'clip_noise_value' must be a finite scalar when clip_noise = 'fixed'")
    }
    if (clip_noise_value < 0) stop("'clip_noise_value' must be non-negative")
    noise_thr <- clip_noise_value
  } else if (clip_noise == "quantile") {
    if (is.null(clip_noise_value) || length(clip_noise_value) != 1 || !is.finite(clip_noise_value)) {
      stop("'clip_noise_value' must be a finite scalar when clip_noise = 'quantile'")
    }
    if (clip_noise_value < 0 || clip_noise_value > 1) {
      stop("'clip_noise_value' must be in [0, 1] when clip_noise = 'quantile'")
    }
    noise_thr <- .abs_quantile(mat_raw, clip_noise_value)
  }

  # Estimate scaling from the RAW assay, not from a clipped intermediate.
  mat_pos_raw <- .split_pos(mat_raw)
  mat_neg_raw <- .split_neg(mat_raw)

  list(
    use_assay  = use_assay,
    transform  = transform,
    noise_thr  = noise_thr,
    weight_pos = .compute_weight(mat_pos_raw, clip_max, clip_max_value),
    weight_neg = .compute_weight(mat_neg_raw, clip_max, clip_max_value)
  )
}

preprocess_assay_apply <- function(lemur_fit_class, spec) {
  mat_raw <- assay(lemur_fit_class, spec$use_assay)

  if (!is.null(spec$noise_thr) && spec$noise_thr > 0) {
    mat_raw <- .apply_fixed_noise_clip(mat_raw, spec$noise_thr)
  }

  mat_pos <- .split_pos(mat_raw)
  mat_neg <- .split_neg(mat_raw)

  if (inherits(mat_pos, "sparseMatrix")) {
    mat_pos <- .apply_transform_sparse(mat_pos, spec$transform, spec$weight_pos)
    mat_neg <- .apply_transform_sparse(mat_neg, spec$transform, spec$weight_neg)
  } else {
    mat_pos <- .apply_transform_dense(mat_pos, spec$transform, spec$weight_pos)
    mat_neg <- .apply_transform_dense(mat_neg, spec$transform, spec$weight_neg)
  }

  mat_pos <- Matrix::drop0(Matrix::Matrix(mat_pos, sparse = TRUE))
  mat_neg <- Matrix::drop0(Matrix::Matrix(mat_neg, sparse = TRUE))

  rn <- rownames(mat_raw)
  if (!is.null(rn)) {
    rn <- make.unique(rn)
    rownames(mat_pos) <- paste0("up|", rn)
    rownames(mat_neg) <- paste0("down|", rn)
  }
  colnames(mat_pos) <- colnames(mat_raw)
  colnames(mat_neg) <- colnames(mat_raw)

  Matrix::rbind2(mat_pos, mat_neg)
}

preprocess_assay <- function(
  lemur_fit_class,
  use_assay,
  transform        = c("none", "tanh", "asinh"),
  clip_noise       = c("none", "fixed", "quantile"),
  clip_noise_value = NULL,
  clip_max         = c("none", "quantile", "fixed"),
  clip_max_value   = NULL
) {
  spec <- preprocess_assay_fit(
    lemur_fit_class = lemur_fit_class,
    use_assay = use_assay,
    transform = transform,
    clip_noise = clip_noise,
    clip_noise_value = clip_noise_value,
    clip_max = clip_max,
    clip_max_value = clip_max_value
  )
  preprocess_assay_apply(lemur_fit_class, spec)
}

preprocess_assay_splits <- function(
  lemur_fit_class,
  use_assay,
  ...
) {
  list(
    train = preprocess_assay(lemur_fit_class$training_data, use_assay, ...),
    test  = preprocess_assay(lemur_fit_class$test_data,     use_assay, ...)
  )
}






compose_multi_condition <- function(split_views, weights = NULL) {

  avail_views <- names(split_views)

  if (is.null(weights)) {
    weights <- rep(1, length(avail_views))
  }
  
  train_views <- lapply(seq_along(split_views), function(i) {
    x <- split_views[[i]]$train * weights[i]
    rownames(x) <- paste0(avail_views[i], "|", rownames(x))
    x
  })

  test_views <- lapply(seq_along(split_views), function(i) {
    x <- split_views[[i]]$test * weights[i]
    rownames(x) <- paste0(avail_views[i], "|", rownames(x))
    x
  })

  list(
    train = do.call(rbind, train_views),
    test  = do.call(rbind, test_views)
  )
}