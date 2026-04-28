#' Calculate Genetic Diversity in a Gridded Bounding Box
#'
#' @param gm A genomaps object containing genetic data and geographic informations, created by [genomaps()] function.
#' @param gmarea Raster file that contains area size of each cell
#' @param bbox Numeric vector of length 4 specifying the bounding box coordinates c(r1, r2, c1, c2)
#' @param revbbox Logical, whether to reverse/invert the bounding box selection. Default FALSE
#'
#' @return A list containing:
#'   \item{N}{Number of samples}
#'   \item{M}{Number of segregating sites}
#'   \item{E}{Number of endemic segregating sites}
#'   \item{thetaw}{Watterson's theta estimate}
#'   \item{thetapi}{Pi (pairwise) diversity estimate}
#'   \item{A}{Total area with data}
#'   \item{Asq}{Area of the bounding box square}
#' @export
#'
#' @examples
#' # Calculate mutation diversity in a 10x10 grid region
#' gmarea <- mar:::.areaofraster(gm1001g$maps$samplemap)
#' div <- mutdiv.gridded(gm1001g, gmarea, bbox = c(1, 10, 2, 8))
mutdiv.gridded <- function(gm, gmarea, bbox, revbbox = FALSE) {
    stopifnot(length(bbox) == 4)
    r1 <- bbox[1]
    r2 <- bbox[2]
    c1 <- bbox[3]
    c2 <- bbox[4]
    nrow <- r2 - r1 + 1
    ncol <- c2 - c1 + 1
    resrow <- terra::res(gm$maps$samplemap)[1]
    rescol <- terra::res(gm$maps$samplemap)[2]
    # calculate area by size of bounding box
    Asq <- .areaofsquare(nrow, ncol, resrow, rescol)
    # if reverse bounding box
    if (revbbox) {
        Asq <- .areaofsquare(dim(gm$maps$samplemap)[1], dim(gm$maps$samplemap)[2], resrow, rescol) - Asq
    }
    # locate cellids from bbox
    cellids <- .rowcol_cellid(gm$maps, bbox, revbbox = revbbox)
    out <- .mutdiv.cellids(gm, gmarea, cellids, Asq)
    return(out)
}

#' Calculate Genetic Diversity for Specified Cells
#'
#' Calculates genetic diversity metrics for a specific set of cells. This function is particularly
#' useful for extinction simulations.
#'
#' @param gm A genomaps object containing genetic data and geographic informations, created by [genomaps()] function.
#' @param gmarea Raster file that contains area size of each cell
#' @param cellids Vector of cell IDs to analyze
#'
#' @return A list containing:
#'   \item{N}{Number of samples}
#'   \item{M}{Number of segregating sites}
#'   \item{E}{Number of endemic segregating sites}
#'   \item{thetaw}{Watterson's theta estimate}
#'   \item{thetapi}{Pi (pairwise) diversity estimate}
#'   \item{A}{Total area with data}
#'   \item{Asq}{Area of the cells being analyzed}
#' @export
#'
#' @examples
#' # Calculate genetic diversity for a specific set of cells
#' gmarea <- mar:::.areaofraster(gm1001g$maps$samplemap)
#' cell_ids <- c(613, 726, 727)
#' div <- mutdiv.cells(gm1001g, gmarea, cellids = cell_ids)
mutdiv.cells <- function(gm, gmarea, cellids) {
    resrow <- terra::res(gm$maps$samplemap)[1]
    rescol <- terra::res(gm$maps$samplemap)[2]
    # calculate area by size of bounding box
    Asq <- .areaofsquare(length(cellids), ncol = 1, resrow, rescol)
    out <- .mutdiv.cellids(gm, gmarea, cellids, Asq)
    return(out)
}

.mutdiv.cellids <- function(gm, gmarea, cellids, Asq) {
    # if no cells
    if (length(cellids) == 0) {
        out <- list(N = NA, M = NA, E = NA, thetaw = NA, thetapi = NA, A = NA, Asq = Asq)
    } else {
        # calculate area by filled raster
        A <- sum(gmarea[cellids])
        # get samples
        sampleids <- .cellid_sample(gm$maps, cellids)
        # calculate genetic diversity
        out <- append(
            .calc_theta(gm, sampleids),
            list(
                A = A,
                Asq = Asq
            )
        )
    }
    return(out)
}

# genetic diversity estimator (use `gm$genotype` matrix)
# ploidy does not matter here. although > diploid is not well-defined.
# TODO: allow L calculations
.calc_theta <- function(gm, sampleid = NULL) {
    ploidy <- .get_genodata(gm$geno, "ploidy")

    if (is.null(sampleid)) {
        sampleid <- 1:(dim(gm$geno$genotype)[2])
    }

    geno <- gm$geno$genotype[, sampleid, drop = FALSE]

    # allele counts
    AC <- matrixStats::rowSums2(geno, na.rm = TRUE)

    # called alleles per site
    n_called_ind <- matrixStats::rowSums2(!is.na(geno))
    called_alleles <- n_called_ind * ploidy

    valid <- called_alleles > 1

    if (!any(valid)) {
        return(list(
            N = length(sampleid), M = 0, E = 0,
            thetaw = 0, thetapi = 0
        ))
    }

    # ---------------------------
    # Theta Pi (pairwise correct)
    # ---------------------------
    num <- sum(
        AC[valid] * (called_alleles[valid] - AC[valid]),
        na.rm = TRUE
    )

    den <- sum(
        choose(called_alleles[valid], 2),
        na.rm = TRUE
    )

    thetapi <- if (den > 0) num / den else 0

    # ---------------------------
    # Theta Watterson
    # ---------------------------
    a_n <- rep(NA_real_, length(called_alleles))
    a_n[valid] <- .Hn(called_alleles[valid] - 1)

    seg_sites <- (AC > 0) & (AC < called_alleles) & valid

    thetaw <- if (any(seg_sites)) {
        sum(seg_sites) / mean(a_n[valid], na.rm = TRUE)
    } else {
        0
    }

    # ---------------------
    # Counts
    # ---------------------
    M <- sum(seg_sites)

    # outgroup handling
    outgeno <- gm$geno$genotype[, -sampleid, drop = FALSE]

    oAC <- matrixStats::rowSums2(outgeno, na.rm = TRUE)
    o_called <- matrixStats::rowSums2(!is.na(outgeno)) * ploidy

    E <- sum(seg_sites & oAC == 0 & o_called > 0)

    return(list(
        N = length(sampleid),
        M = M,
        E = E,
        thetaw = thetaw,
        thetapi = thetapi
    ))
}
