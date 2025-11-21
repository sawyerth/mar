test_that("MARsampling works with gm1001g genomaps data", {
    data("gm1001g", package = "mar")

    result <- mar::MARsampling(
        gm = gm1001g,
        scheme = "random",
        nrep = 2,
        quorum = FALSE,
        myseed = 123
    )

    expect_s3_class(result, "marsamp")
    expect_gt(nrow(result), 0)
    expect_identical(attr(result, "scheme"), "random")
    expect_true(all(c("N", "M", "E", "thetaw", "thetapi", "A", "Asq", "extent") %in% names(result)))
    expect_false(any(is.na(result$Asq)))
})
