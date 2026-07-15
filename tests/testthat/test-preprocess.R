test_that("preprocess_assay_fit returns valid specification", {

    fit <- tiny_lemur_fit()

    spec <- preprocess_assay_fit(
        fit,
        use_assay = "DE"
    )

    expect_type(spec, "list")

    expect_named(
        spec,
        c(
            "use_assay",
            "transform",
            "noise_thr",
            "weight_pos",
            "weight_neg"
        )
    )

    expect_equal(spec$use_assay, "DE")

})


test_that("invalid assay throws informative error", {

    fit <- tiny_lemur_fit()

    expect_error(

        preprocess_assay_fit(
            fit,
            use_assay = "foobar"
        ),

        "use_assay"

    )

})


test_that("quantile clipping computes threshold", {

    fit <- tiny_lemur_fit()

    spec <- preprocess_assay_fit(

        fit,

        use_assay = "DE",

        clip_noise = "quantile",

        clip_noise_value = 0.2

    )

    expect_true(
        spec$noise_thr >= 0
    )

})





test_that("preprocess_assay_apply doubles features", {

    fit <- tiny_lemur_fit()

    spec <- preprocess_assay_fit(
        fit,
        "DE"
    )

    mat <- preprocess_assay_apply(
        fit,
        spec
    )

    expect_equal(

        nrow(mat),

        2 * nrow(fit)

    )

    expect_equal(

        ncol(mat),

        ncol(fit)

    )

})




test_that("rows receive up/down prefixes", {

    fit <- tiny_lemur_fit()

    spec <- preprocess_assay_fit(
        fit,
        "DE"
    )

    mat <- preprocess_assay_apply(
        fit,
        spec
    )

    expect_true(
        all(
            startsWith(
                rownames(mat)[1:nrow(fit)],
                "up|"
            )
        )
    )

})





test_that("output is sparse", {

    fit <- tiny_lemur_fit()

    spec <- preprocess_assay_fit(
        fit,
        "DE"
    )

    mat <- preprocess_assay_apply(
        fit,
        spec
    )

    expect_s4_class(
        mat,
        "dgCMatrix"
    )

})





test_that("preprocess_assay returns train/test", {

    fit <- tiny_lemur_fit()

    out <- preprocess_assay(
        fit,
        "DE"
    )

    expect_named(
        out,
        c("train","test")
    )

    expect_s4_class(
        out$train,
        "dgCMatrix"
    )

    expect_s4_class(
        out$test,
        "dgCMatrix"
    )

})





test_that("compose_multi_view combines views", {

    fit <- tiny_lemur_fit()

    de <- preprocess_assay(
        fit,
        "DE"
    )

    counts <- preprocess_assay(
        fit,
        "counts"
    )

    out <- compose_multi_view(
        list(
            DE = de,
            counts = counts
        )
    )

    expect_named(
        out,
        c("train","test")
    )

})






test_that("compose_multi_view prefixes views", {

    fit <- tiny_lemur_fit()

    de <- preprocess_assay(fit, "DE")
    counts <- preprocess_assay(fit, "counts")

    out <- compose_multi_view(
        list(
            DE = de,
            counts = counts
        )
    )

    expect_true(
        any(startsWith(rownames(out$train), "DE|"))
    )

    expect_true(
        any(startsWith(rownames(out$train), "counts|"))
    )

})