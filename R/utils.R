.parse_feature_names <- function(x) {

  parts <- strsplit(x, "|", fixed = TRUE)

  info <- lapply(seq_along(parts), function(i) {

    p <- parts[[i]]

    if (length(p) >= 3 && p[2] %in% c("up", "down")) {

      data.frame(
        id = x[i],
        view = p[1],
        direction = p[2],
        feature = paste(p[-c(1, 2)], collapse = "|"),
        stringsAsFactors = FALSE
      )

    } else if (length(p) >= 2 && p[1] %in% c("up", "down")) {

      data.frame(
        id = x[i],
        view = NA_character_,
        direction = p[1],
        feature = paste(p[-1], collapse = "|"),
        stringsAsFactors = FALSE
      )

    } else {

      data.frame(
        id = x[i],
        view = NA_character_,
        direction = NA_character_,
        feature = paste(p, collapse = "|"),
        stringsAsFactors = FALSE
      )

    }
  })

  do.call(rbind, info)
}
