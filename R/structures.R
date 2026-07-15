.init_group <- function() {
  structure(
    list(
      elements = character(0),
      loadings = numeric(0),
      metadata = list()
    ),
    class = "grp"
  )
}

.init_factor <- function() {
  structure(
    list(
      views = list(),
      cells = .init_group(),
      test  = NULL,
      metadata = list()
    ),
    class = "biclust_factor"
  )
}


.init_BiclustResult <- function() {
  structure(
    list(
      consensus       = list(),
      factors         = list(),
      biclusters      = list(),
      summary         = list(),
      diagnostics     = list(),
      analyses        = list(),
      metadata        = list()
    ),
    class = "BiclustResult"
  )
}




.init_consensus <- function() {
    structure(
        list(
            models = list(),
            alignments = NULL,

            gene_support = NULL,
            cell_support = NULL,

            gene_loading_mean = NULL,
            cell_loading_mean = NULL,

            sample_names = NULL,
            feature_names = NULL,

            metadata = list()
        ),
        class = "biclust_consensus"
    )
}