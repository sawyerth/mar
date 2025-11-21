test_that(".valid_lonlat enforces matrix structure", {
    good <- matrix(c(1, 2, 3, 4), ncol = 2)
    expect_invisible(mar:::.valid_lonlat(good))

    expect_error(mar:::.valid_lonlat(list(1, 2)))
    expect_error(mar:::.valid_lonlat(matrix(1:3, ncol = 3)))
    expect_error(mar:::.valid_lonlat(matrix(c(1, NA), ncol = 2)))
})

test_that(".areaofsquare multiplies dimensions", {
    expect_equal(mar:::.areaofsquare(2, 3, 4, 5), 120)
})

test_that(".rowcol_extent creates the correct extent", {
    gm <- gm1001g
    bbox <- c(1, 10, 1, 10)
    expected <- raster::extent(gm$maps$samplemap, bbox[1], bbox[2], bbox[3], bbox[4])
    result <- mar:::.rowcol_extent(gm$maps, bbox)
    expect_equal(result, expected)
})
