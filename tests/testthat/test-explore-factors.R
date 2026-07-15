test_that("explore_factors returns factor_exploration", {

    mat <- tiny_preprocessed()

    out <- explore_factors(
        mat,
        ks = 2:3,
        n_restarts = 2
    )

    expect_s3_class(
        out,
        "factor_exploration"
    )

})



test_that("explore_factors returns factor_exploration", {

    mat <- tiny_preprocessed()

    out <- explore_factors(
        mat,
        ks = 2:3,
        n_restarts = 2
    )

    expect_s3_class(
        out,
        "factor_exploration"
    )

})


test_that("invalid k throws error", {

    mat <- tiny_preprocessed()

    expect_error(

        explore_factors(
            mat,
            ks = 0
        ),

        "positive integers"

    )

})


test_that("invalid test fraction throws error", {

    mat <- tiny_preprocessed()

    expect_error(

        explore_factors(
            mat,
            ks = 2:3,
            test_fraction = 2
        ),

        "test_fraction"

    )

})


test_that("print works", {

    mat <- tiny_preprocessed()

    out <- explore_factors(
        mat,
        ks = 2:3,
        n_restarts = 1
    )
  
  expect_output(
      print(out)
  )

})



test_that("plot returns invisibly", {

    mat <- tiny_preprocessed()

    out <- explore_factors(
        mat,
        ks = 2:3,
        n_restarts = 1
    )

    pdf(NULL)

    on.exit(dev.off())

    expect_output(
        print(out)
    )

})

