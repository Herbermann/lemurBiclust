.TRANSFORMS <- list(
  tanh  = function(w = 1, ...) function(x) tanh(w * x) / w,
  none    = function(...)  identity,
  clip    = function(threshold) function(x) pmax(x, threshold)
)

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

  # --- dispatch ---

  mat_raw <- assay(lemur_fit_class, use_assay)
  is_sparse <- is(mat_raw, "sparseMatrix")

  # --- compute weight (saturation scale) ---
  weight <- if (!is.null(clip_max_value)) {
    if (clip_max == "quantile") {
      1 / quantile(mat_raw, clip_max_value, na.rm = TRUE)
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

  # --- compute threshold (noise floor) ---
  thrsh <- if (!is.null(clip_noise_value)) {
    if (clip_noise == "fixed") {
      clip_noise_value
    } else {
      0
    }
  } else {
    if (clip_noise != "none")
      warning("no valid clip_noise_value provided, defaulting to thrsh = 0")
    0
  }

  # --- apply noise floor clipping ---
  # pmax introduces new zeros, so drop them before transform for efficiency
  if (thrsh != 0) {
    mat_raw <- .TRANSFORMS[["clip"]](threshold = thrsh)(mat_raw)
    if (is_sparse)
      mat_raw <- Matrix::drop0(mat_raw)
  }

  # --- apply transform ---
  transform_fn <- .TRANSFORMS[[transform]](w = weight)
  mat <- transform_fn(mat_raw)

  # --- restore sparsity if input was sparse ---
  if (is_sparse && !is(mat, "sparseMatrix"))
    mat <- Matrix::Matrix(mat, sparse = TRUE)

  # --- split into positive and negative parts ---
  mat_pos <- pmax(mat,  0)
  mat_neg <- pmax(-mat, 0)

  # ensure sparse after split — pmax densifies, and we want zeros dropped
  mat_pos <- Matrix::drop0(Matrix::Matrix(mat_pos, sparse = TRUE))
  mat_neg <- Matrix::drop0(Matrix::Matrix(mat_neg, sparse = TRUE))

  rownames(mat_pos) <- paste0("up_",   make.unique(rownames(mat)))
  rownames(mat_neg) <- paste0("down_", make.unique(rownames(mat)))
  colnames(mat_pos) <- colnames(mat)
  colnames(mat_neg) <- colnames(mat)

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