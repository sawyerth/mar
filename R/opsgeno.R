# matrix validators. Also validates that there are no invariant sites
.valid_genotype <- function(genotype, ploidy) {
    if (!is.matrix(genotype)) {
        stop("genotype must be a matrix")
    }

    valid_vars <- seq(0, ploidy)

    all_vars <- unique(as.vector(genotype))

    # allow NA, but reject invalid allele states
    badvars <- all_vars[!is.na(all_vars) & !(all_vars %in% valid_vars)]

    if (length(badvars) > 0) {
        stop(paste0(
            "genotype contains ",
            length(badvars),
            " unique invalid values: ",
            toString(head(badvars)), "..."
        ))
    }

    return(invisible())
}

# genotype operations that works on both margeno and data frame objects
.get_genodata <- function(x, what = c("sample.id", "variant.id", "position", "chromosome", "genotype", "ploidy", "num.variant")) {
    what <- match.arg(what)
    if ("margeno" %in% class(x)) {
        return(switch(what,
            "num.variant" = length(x$variant.id),
            x[[what]]
        ))
    } else if (is.data.frame(x)) {
        return(switch(what,
            "sample.id"   = names(x)[-(1:9)],
            "variant.id"  = x$ID,
            "chromosome"  = x$CHROM,
            "position"    = x$POS,
            "ploidy"      = .infer_ploidy(x),
            "num.variant" = nrow(x),
            "genotype"    = .vcfgt2mat(x)
        ))
    } else {
        stop("x must be either a data frame (VCF) or a margeno object")
    }
}

# Extract genotype dosage matrix (variant x sample) from a VCF data frame
.vcfgt2mat <- function(x) {
    sample_cols <- x[, -(1:9), drop = FALSE]
    gt <- lapply(sample_cols, function(col) {
        # Extract GT field (before first colon)
        gt_field <- sub(":.*", "", col)
        vapply(gt_field, function(g) {
            if (g %in% c(".", "./.", ".|.")) {
                return(NA_integer_)
            }
            alleles <- as.integer(strsplit(g, "[/|]")[[1]])
            sum(alleles)
        }, integer(1))
    })
    matrix(unlist(gt),
        nrow = nrow(x), ncol = length(gt),
        dimnames = list(x$ID, names(gt))
    )
}

.infer_ploidy <- function(x) {
    sample_cols <- x[, -(1:9), drop = FALSE]
    first_gt <- sub(":.*", "", sample_cols[[1]][1])
    if (first_gt %in% c(".", "./.", ".|.")) {
        return(2L)
    } # fallback for missing
    length(strsplit(first_gt, "[/|]")[[1]])
}

# convert seqGetData output to genotype matrix in margeno
# in SeqArray, dim1 = ploidy, dim2 = number of sample, dim3 = number of variant
# in margeno$genotype, dim1 = number of variant (row), dim2 = number of sample (col), only diploid or pseudo-haploid allowed
.seqgeno2mat <- function(gg, ploidy) {
    stopifnot(class(gg) == "array")
    stopifnot(length(dim(gg)) == 3 & dim(gg)[1] == ploidy)
    # sum up the allele counts
    df <- apply(gg, c(3, 2), sum)
    attr(df, "dimnames") <- NULL
    .valid_genotype(df, ploidy)
    return(df)
}

# convert data frame object to margeno object
.dataframe2margeno <- function(x) {
    stopifnot(is.data.frame(x))

    # create margeno object
    out <- .new_margeno(
        sample.id = .get_genodata(x, "sample.id"),
        variant.id = .get_genodata(x, "variant.id"),
        position = .get_genodata(x, "position"),
        chromosome = .get_genodata(x, "chromosome"),
        genotype = .get_genodata(x, "genotype"),
        ploidy = .get_genodata(x, "ploidy")
    )
    return(out)
}

# # convert genotype matrix to SeqArray format (3D array)
# .mat2seqgeno <- function(gg) {
#     stopifnot(class(gg) == "matrix")
#     stopifnot(length(dim(gg)) == 2)
#     df <- array(dim=c(2, rev(dim(gg))))
#     df[1,,] <- apply(gg, 1, function(xx) as.integer(xx > 1))
#     df[2,,] <- apply(gg, 1, function(xx) as.integer(xx > 0))
#     return(df)
# }

# # convert margeno object to SeqArray object (TODO: debug)
# # inspired by: https://github.com/zhengxwen/SeqArray/issues/62
# .margeno2seqarray <- function(x, gds.fn = NULL, opengds = FALSE) {
#     stopifnot(class(x) == "margeno")

#     # Create a new GDS file
#     if(is.null(gds.fn)) {
#         gds.fn <- tempfile(fileext = ".gds")
#     }

#     gds <- gdsfmt::createfn.gds(gds.fn)

#     # Add the data
#     gdsfmt::add.gdsn(gds, "sample.id", x$sample.id)
#     gdsfmt::add.gdsn(gds, "variant.id", x$variant.id)
#     gdsfmt::add.gdsn(gds, "position", x$position)
#     gdsfmt::add.gdsn(gds, "chromosome", x$chromosome)
#     gdsfmt::add.gdsn(gds, "allele", rep("N,N", length(x$variant.id)))
#     geno_array <- .mat2seqgeno(x$genotype)
#     gdsfmt::add.gdsn(gds, "genotype", geno_array)

#     gdsfmt::put.attr.gdsn(gds$root, "FileFormat", "SEQ_ARRAY")
#     gdsfmt::put.attr.gdsn(gds$root, "FileVersion", "v1.0")
#     gdsfmt::addfolder.gdsn(gds, "description")

#     # Close the GDS file
#     gdsfmt::closefn.gds(gds)

#     # Open as SeqVarGDSClass
#     if (opengds) {
#         seqfile <- SeqArray::seqOpen(gds.fn)
#         return(seqfile)
#     } else {
#         return(gds.fn)
#     }
# }

.filter_genosample <- function(x, sample.id) {
    stopifnot(all(sample.id %in% .get_genodata(x, "sample.id")))
    if ("margeno" %in% class(x)) {
        idx <- match(sample.id, x$sample.id)
        x$sample.id <- x$sample.id[idx]
        x$genotype <- as.matrix(x$genotype[, idx])
        return(x)
    } else if (is.data.frame(x)) {
        sample_cols <- names(x)[-(1:9)]
        idx <- match(sample.id, sample_cols)
        return(x[, c(1:9, 9 + idx)])
    } else {
        stop("x must be either a data frame (VCF) or a margeno object")
    }
}

.filter_genovariant <- function(x, variant.id) {
    stopifnot(all(variant.id %in% .get_genodata(x, "variant.id")))
    if ("margeno" %in% class(x)) {
        idx <- match(variant.id, x$variant.id)
        x$variant.id <- x$variant.id[idx]
        x$position <- x$position[idx]
        x$chromosome <- x$chromosome[idx]
        x$genotype <- as.matrix(x$genotype[idx, ])
        return(x)
    } else if (is.data.frame(x)) {
        idx <- match(variant.id, x$ID)
        return(x[idx, ])
    } else {
        stop("x must be either a data frame (VCF) or a margeno object")
    }
}

# author Feng Li, Department of Statistics, Stockholm University, Sweden.
# https://github.com/feng-li/flutils/blob/master/R/math__hamonic.R
.Hn <- function(n) {
    out <- digamma(n + 1) - digamma(1)
    return(out)
}
