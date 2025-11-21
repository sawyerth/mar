cleanup_workdir <- function(path) {
    if (dir.exists(path)) {
        unlink(path, recursive = TRUE, force = TRUE)
    }
}

test_that("MARPIPELINE processes example extdata", {
    name <- "extdata_example"
    workdir <- tempfile("marpipeline_extdata_")
    dir.create(workdir)
    on.exit(cleanup_workdir(workdir), add = TRUE)

    genofile <- system.file("extdata", "genome.tsv", package = "mar")
    lonlatfile <- system.file("extdata", "lonlat.csv", package = "mar")

    oldwd <- getwd()
    on.exit(setwd(oldwd), add = TRUE)

    expect_invisible(MARPIPELINE(
        name = name,
        workdir = workdir,
        genofile = genofile,
        lonlatfile = lonlatfile,
        filetype = "text"
    ))
})
