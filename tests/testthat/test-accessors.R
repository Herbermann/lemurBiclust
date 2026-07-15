test_that("biclusters returns a BiclusterList", {

    fit <- tiny_biclust_result()

    bcs <- biclusters(fit)

    expect_s3_class(
        bcs,
        "BiclusterList"
    )

    expect_length(
        bcs,
        2
    )

})


test_that("single bicluster can be extracted by index", {

    fit <- tiny_biclust_result()

    bc <- biclusters(fit, 1)

    expect_true(
        is.list(bc)
    )

    expect_true(
        "cells" %in% names(bc)
    )

})



test_that("single bicluster can be extracted by id", {

    fit <- tiny_biclust_result()

    id <- names(biclusters(fit))[1]

    bc <- biclusters(
        fit,
        bic = id
    )

    expect_true(
        "cells" %in% names(bc)
    )

})



test_that("invalid bic id throws error", {

    fit <- tiny_biclust_result()

    expect_error(

        biclusters(
            fit,
            "foobar"
        ),

        "Unknown bicluster"

    )

})



test_that("invalid object throws error", {

    expect_error(

        biclusters(1),

        "BiclustResult"

    )

})



test_that("factors accessor works", {

    fit <- tiny_biclust_result()

    expect_length(
        factors(fit),
        2
    )

})



test_that("consensus accessor works", {

    fit <- tiny_biclust_result()

    expect_true(
        is.list(
            consensus(fit)
        )
    )

})



test_that("analyses is initially empty", {

    fit <- tiny_biclust_result()

    expect_length(
        analyses(fit),
        0
    )

})


test_that("diagnostics accessor works", {

    fit <- tiny_biclust_result()

    expect_true(
        is.list(
            diagnostics(fit)
        )
    )

})



test_that("composition works for one label", {

    fit <- tiny_biclust_result()

    comp <- bicluster_composition(

        biclusters(fit),

        tiny_lemur_fit(),

        label = "cell_type"

    )

    expect_true(
        all(
            c(
                "bicluster",
                "variable",
                "label",
                "n",
                "fraction_bic",
                "fraction_label"
            ) %in%
            colnames(comp)
        )
    )

})



test_that("composition supports multiple labels", {

    fit <- tiny_biclust_result()

    comp <- bicluster_composition(

        biclusters(fit),

        tiny_lemur_fit(),

        label = c(
            "cell_type",
            "condition"
        )

    )

    expect_equal(

        sort(unique(comp$variable)),

        c(
            "cell_type",
            "condition"
        )

    )

})




test_that("composition supports both train and test", {

    fit <- tiny_biclust_result()

    comp <- bicluster_composition(

        biclusters(fit),

        tiny_lemur_fit(),

        label = "cell_type",

        cell_slot = "both"

    )

    expect_equal(

        sort(unique(comp$cell_slot)),

        c(
            "cells",
            "test"
        )

    )

})



test_that("composition supports single bicluster", {

    fit <- tiny_biclust_result()

    bc <- biclusters(
        fit,
        1
    )

    comp <- bicluster_composition(

        bc,

        tiny_lemur_fit(),

        label = "cell_type"

    )

    expect_true(
        "bicluster" %in%
        names(comp)
    )

})