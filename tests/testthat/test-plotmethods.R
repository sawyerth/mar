context("plotmethods")

test_that("plot methods work with gm1001g", {
    data(gm1001g)

    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off())

    expect_invisible(plot(gm1001g$maps))

    sfs_obj <- sfs(
        AC = mar:::.get_AC(gm1001g$geno),
        N = length(gm1001g$geno$sample.id),
        ploidy = gm1001g$geno$ploidy
    )
    expect_invisible(plot(sfs_obj))

    marsamp_df <- MARsampling(gm = gm1001g, scheme = "random", nrep = 1, xfrac = 0.1, quorum = FALSE, myseed = 7)
    expect_s3_class(marsamp_df, "marsamp")
    expect_invisible(plot(marsamp_df))

    extinction_df <- MARextinction(gm = gm1001g, scheme = "random", nrep = 1, xfrac = 0.1, myseed = 7)
    expect_s3_class(extinction_df, "marextinct")
    expect_invisible(plot(extinction_df))
})
