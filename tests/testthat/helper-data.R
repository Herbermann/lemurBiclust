progressr::handlers("void")

tiny_counts <- function() {

    set.seed(1)

    matrix(
        rpois(20 * 120, lambda = 5),
        nrow = 20,
        dimnames = list(
            paste0("gene", seq_len(20)),
            paste0("cell", seq_len(120))
        )
    )
}


tiny_sce <- function() {

    counts <- tiny_counts()

    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(
            counts = counts,
            logcounts = log1p(counts)
        )
    )

    SummarizedExperiment::colData(sce)$sample <-
        rep(paste0("S", 1:6), each = 20)

    SummarizedExperiment::colData(sce)$condition <-
        rep(c("A", "B"), each = 60)

    SummarizedExperiment::colData(sce)$individual <-
        rep(paste0("I", 1:12), each = 10)

    SummarizedExperiment::colData(sce)$cell_type <-
        sample(
            c("T", "B", "NK", "Mono"),
            ncol(sce),
            replace = TRUE
        )

    sce

}


tiny_lemur_fit <- local({

    fit <- NULL

    function() {

        if (!is.null(fit))
            return(fit)

        sce <- tiny_sce()
        set.seed(1)
        fit <- lemur::lemur(
            sce,
            design = ~ condition,
            test_fraction = 0.3,
            n_embedding = 4,
            verbose = FALSE
        )
        fit <- lemur::test_de(fit, contrast = cond(condition = "A") -  cond(condition = "B") )
        fit

    }
})



tiny_preprocessed <- function(
    assay = "DE",
    ...
) {

    preprocess_assay(
        tiny_lemur_fit(),
        use_assay = assay,
        ...
    )

}



tiny_biclust_result <- local({

    fit <- NULL

    function() {

        if (!is.null(fit))
            return(fit)

        fit <<- biclust(
            tiny_preprocessed(),
            k = 2,
            n_reps = 2,
            threads = 1,
            threads_RcppML = 1,
            verbose = FALSE
        )

        fit

    }

})