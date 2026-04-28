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

test_that("MARPIPELINE runs with missing genotype data (extdata)", {
    name <- "extdata_missing"
    workdir <- tempfile("marpipeline_missing_")
    dir.create(workdir)
    on.exit(cleanup_workdir(workdir), add = TRUE)

    genofile <- system.file("extdata", "genome.tsv", package = "mar")
    lonlatfile <- system.file("extdata", "lonlat.csv", package = "mar")

    # -----------------------------
    # Load and inject missing data
    # -----------------------------
    geno <- read.table(genofile, header = TRUE, stringsAsFactors = FALSE)

    geno <- data.matrix(geno)

    # ONLY sample valid numeric entries
    valid_idx <- which(geno %in% c(0, 1, 2))

    set.seed(123)
    geno[sample(valid_idx, 5)] <- NA

    temp_genofile <- file.path(workdir, "genome_na.tsv")
    write.table(geno,
        temp_genofile,
        sep = "\t",
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE,
        na = "NA"
    )


    oldwd <- getwd()
    on.exit(setwd(oldwd), add = TRUE)

    # -----------------------------
    # Run pipeline (minimal steps)
    # -----------------------------
    expect_invisible(
        MARPIPELINE(
            name = name,
            workdir = workdir,
            genofile = temp_genofile,
            lonlatfile = lonlatfile,
            filetype = "text",
            option_geno = list(ploidy = 2, maxsnps = 100000),
            option_map = list(mapres = NULL, mapcrs = "OGC:CRS84"),
            option_sadsfs = list(sad_models = .sad_models, folded = TRUE),
            option_marext = list(
                scheme = "random",
                nrep = 1,
                xfrac = 0.5,
                quorum = FALSE,
                animate = FALSE
            ),
            marsteps = c("data", "gm"), # keep lightweight
            saveobj = FALSE
        )
    )
})
