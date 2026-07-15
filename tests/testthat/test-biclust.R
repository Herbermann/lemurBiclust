test_that("biclust returns a BiclustResult", {

    mat <- tiny_preprocessed()

    fit <- biclust(
        mat,
        k = 2,
        n_reps = 2,
        verbose = FALSE
    )

    expect_s3_class(
        fit,
        "BiclustResult"
    )

})


test_that("biclust populates result slots", {

    mat <- tiny_preprocessed()

    fit <- biclust(
        mat,
        k = 2,
        n_reps = 2,
        verbose = FALSE
    )

    expect_true(length(fit$consensus) > 0)
    expect_true(length(fit$factors) == 2)
    expect_true(length(fit$biclusters) == 2)

})



test_that("one bicluster per factor", {

    mat <- tiny_preprocessed()

    fit <- biclust(
        mat,
        k = 2,
        n_reps = 2,
        verbose = FALSE
    )

    expect_equal(
        length(factors(fit)),
        length(biclusters(fit))
    )

})



test_that("test projection can be disabled", {

    mat <- tiny_preprocessed()

    fit <- biclust(
        mat,
        k = 2,
        n_reps = 2,
        use_test = FALSE,
        verbose = FALSE
    )

    bcs <- biclusters(fit)

    expect_false(
        "test" %in% names(bcs[[1]])
    )

})


test_that("test projection is available", {

    mat <- tiny_preprocessed()

    fit <- biclust(
        mat,
        k = 2,
        n_reps = 2,
        use_test = TRUE,
        verbose = FALSE
    )

    bcs <- biclusters(fit)

    expect_true(
        all(
            vapply(
                bcs,
                function(x) "test" %in% names(x),
                logical(1)
            )
        )
    )

    expect_true(
        all(
            vapply(
                bcs,
                function(x) length(x$test) > 0,
                logical(1)
            )
        )
    )

})


test_that("parameters are recorded", {

    mat <- tiny_preprocessed()

    fit <- biclust(
        mat,
        k = 2,
        n_reps = 2,
        verbose = FALSE
    )

    expect_equal(
        fit$metadata$gene_support_threshold,
        0.9
    )

    expect_equal(
        fit$metadata$cell_support_threshold,
        0.9
    )

})



test_that("biclusters contain unique genes", {

    mat <- tiny_preprocessed()

    fit <- biclust(
        mat,
        k = 2,
        n_reps = 2,
        verbose = FALSE
    )

    bcs <- biclusters(fit)

    for (bc in bcs) {

        views <- setdiff(
            names(bc),
            c("cells", "test")
        )

        for (v in views) {

        expect_identical(
            anyDuplicated(names(bc[[v]])),
            0L
        )

        }

    }

})



test_that("changing support threshold changes biclusters", {

    mat <- tiny_preprocessed()

    fit1 <- biclust(
        mat,
        k = 2,
        gene_support_thrs = 0.5,
        n_reps = 2,
        verbose = FALSE
    )

    fit2 <- biclust(
        mat,
        k = 2,
        gene_support_thrs = 0.95,
        n_reps = 2,
        verbose = FALSE
    )

    expect_true(
        length(biclusters(fit1)[[1]]$default) >=
        length(biclusters(fit2)[[1]]$default)
    )

})




test_that("complete biclustering pipeline works", {

    fit <- biclust(
        tiny_preprocessed(),
        k = 2,
        n_reps = 2,
        verbose = FALSE
    )

    expect_s3_class(
        fit,
        "BiclustResult"
    )

    expect_length(
        biclusters(fit),
        2
    )

    expect_length(
        factors(fit),
        2
    )

    expect_length(
        analyses(fit),
        0
    )

})




test_that("biclust supports multiple views", {

    fit <- tiny_lemur_fit()

    mat <- compose_multi_view(
        list(
            "view1" = preprocess_assay(
                fit,
                use_assay = "DE"
            ),
            "view2" = preprocess_assay(
                fit,
                use_assay = "DE"
            )
        )
    )

    res <- biclust(
        mat,
        k = 2,
        n_reps = 2
    )

    bcs <- biclusters(res)

    expect_length(
        bcs,
        2
    )

    expect_true(
        all(
            vapply(
                bcs,
                function(x)
                    all(c("view1", "view2") %in%
                            setdiff(names(x), c("cells", "test"))),
                logical(1)
            )
        )
    )

})