.TRANSFORMS <- list(
  tanh  = function(w = 1, ...) function(x) tanh(w * x) / w,
  none    = function(...)  identity,
  clip    = function(threshold) function(x) pmax(x, threshold)
)

.apply_fixed_noise_clip <- function(mat, thrsh) {
  if (thrsh <= 0) return(mat)

  if (inherits(mat, "sparseMatrix")) {
    if (length(mat@x)) {
      mat@x[abs(mat@x) < thrsh] <- 0
    }
    Matrix::drop0(mat)
  } else {
    mat[abs(mat) < thrsh] <- 0
    mat
  }
}

preprocess_assay <- function(
  lemur_fit_class,
  use_assay,
  transform        = c("none", "tanh", "asinh"),
  clip_noise       = c("none", "fixed"),
  clip_noise_value = NULL,
  clip_max         = c("none", "quantile", "fixed"),
  clip_max_value   = NULL
) {
  # --- input validation ---
  transform  <- match.arg(transform)
  clip_noise <- match.arg(clip_noise)
  clip_max   <- match.arg(clip_max)

  if (!use_assay %in% assayNames(lemur_fit_class))
    stop(sprintf(
      "'use_assay' must be one of: %s",
      paste(assayNames(lemur_fit_class), collapse = ", ")
    ))

  mat_raw   <- assay(lemur_fit_class, use_assay)
  is_sparse <- is(mat_raw, "sparseMatrix")

  # --- apply noise floor clipping ---
  if (clip_noise == "fixed") {
    if (is.null(clip_noise_value) || length(clip_noise_value) != 1 || !is.finite(clip_noise_value)) {
      stop("'clip_noise_value' must be a finite scalar when clip_noise = 'fixed'")
    }
    thrsh <- clip_noise_value
    if (thrsh < 0) stop("'clip_noise_value' must be non-negative")
    mat_raw <- .apply_fixed_noise_clip(mat_raw, thrsh)
  }

  # --- split into positive and negative parts first ---
  mat_pos_raw <- pmax(mat_raw,  0)
  mat_neg_raw <- pmax(-mat_raw, 0)

  if (is_sparse) {
    mat_pos_raw <- Matrix::drop0(Matrix::Matrix(mat_pos_raw, sparse = TRUE))
    mat_neg_raw <- Matrix::drop0(Matrix::Matrix(mat_neg_raw, sparse = TRUE))
  }

  # --- compute separate weights per half ---
  .compute_weight <- function(mat_half, clip_max, clip_max_value) {
    if (!is.null(clip_max_value)) {
      if (clip_max == "quantile") {
        nz <- mat_half[mat_half > 0]
        if (length(nz) == 0) return(1)
        1 / quantile(nz, clip_max_value, na.rm = TRUE)
      } else if (clip_max == "fixed") {
        1 / clip_max_value
      } else {
        1
      }
    } else {
      if (clip_max != "none")
        warning("no valid clip_max_value provided, defaulting to weight = 1")
      1
    }
  }

  weight_pos <- .compute_weight(mat_pos_raw, clip_max, clip_max_value)
  weight_neg <- .compute_weight(mat_neg_raw, clip_max, clip_max_value)

  # --- apply transform to each half independently ---
  transform_fn_pos <- .TRANSFORMS[[transform]](w = weight_pos)
  transform_fn_neg <- .TRANSFORMS[[transform]](w = weight_neg)

  mat_pos <- transform_fn_pos(mat_pos_raw)
  mat_neg <- transform_fn_neg(mat_neg_raw)

  # --- ensure sparse ---
  mat_pos <- Matrix::drop0(Matrix::Matrix(mat_pos, sparse = TRUE))
  mat_neg <- Matrix::drop0(Matrix::Matrix(mat_neg, sparse = TRUE))

  # --- set dimnames ---
  rownames(mat_pos) <- paste0("up_",   make.unique(rownames(mat_raw)))
  rownames(mat_neg) <- paste0("down_", make.unique(rownames(mat_raw)))
  colnames(mat_pos) <- colnames(mat_raw)
  colnames(mat_neg) <- colnames(mat_raw)

  Matrix::rbind2(mat_pos, mat_neg)
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