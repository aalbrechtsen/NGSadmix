#!/usr/bin/env Rscript

usage <- function() {
    cat(
        "Usage:\n",
        "  Rscript scripts/Qconv.R <likes_file> <qlist_file> <threshold> [--count-only]\n",
        sep = ""
    )
}

fail <- function(message, status = 1) {
    cat(sprintf("[ERROR] %s\n", message), file = stderr())
    quit(save = "no", status = status)
}

args <- commandArgs(trailingOnly = TRUE)

if (!(length(args) %in% c(3, 4))) {
    usage()
    quit(save = "no", status = 2)
}

likes_file <- args[1]
qlist_file <- args[2]
threshold <- suppressWarnings(as.numeric(args[3]))
count_only <- length(args) == 4 && args[4] == "--count-only"

if (!file.exists(likes_file)) {
    fail(sprintf("likes file not found: %s", likes_file), 3)
}
if (!file.exists(qlist_file)) {
    fail(sprintf("Q-list file not found: %s", qlist_file), 3)
}
if (is.na(threshold) || threshold < 0) {
    fail("threshold must be a non-negative number", 4)
}

read_vector_file <- function(path) {
    dat <- tryCatch(read.table(path, header = FALSE, stringsAsFactors = FALSE),
                    error = function(e) fail(sprintf("could not read %s: %s", path, e$message), 5))
    if (nrow(dat) == 0) {
        fail(sprintf("empty file: %s", path), 5)
    }
    dat
}

extract_ll_vector <- function(path) {
    dat <- read_vector_file(path)
    if (ncol(dat) == 1) {
        values <- dat[[1]]
    } else {
        values <- dat[[ncol(dat)]]
    }
    values
}

read_q_matrix <- function(path) {
    if (!file.exists(path)) {
        fail(sprintf("Q matrix file not found: %s", path), 6)
    }
    mat <- tryCatch(as.matrix(read.table(path, header = FALSE)),
                    error = function(e) fail(sprintf("could not read Q matrix %s: %s", path, e$message), 6))
    if (length(dim(mat)) != 2 || nrow(mat) == 0 || ncol(mat) == 0) {
        fail(sprintf("invalid Q matrix shape in %s", path), 6)
    }
    t(mat)
}

get_fast <- function(Q, Qold) {
    npop <- nrow(Qold)
    if (nrow(Q) != npop || ncol(Q) != ncol(Qold)) {
        fail("Q matrices do not have matching dimensions", 7)
    }

    res <- matrix(NA_real_, nrow = nrow(Q), ncol = 2)
    for (g in seq_len(nrow(Q))) {
        w <- rowSums((matrix(Q[g, ], nrow = npop, ncol = ncol(Qold), byrow = TRUE) - Qold)^2)
        res[g, ] <- c(which.min(w), min(w))
    }

    dup_targets <- unique(res[duplicated(res[, 1]), 1])
    for (target in dup_targets) {
        idx <- which(res[, 1] == target)
        worst <- idx[which.max(res[idx, 2])]
        res[worst, 1] <- npop + 1
    }

    Q[order(res[, 1]), , drop = FALSE]
}

ll <- suppressWarnings(as.numeric(extract_ll_vector(likes_file)))
if (any(is.na(ll))) {
    fail("likes file contains non-numeric values", 5)
}

qlist <- read_vector_file(qlist_file)[[1]]
if (length(ll) != length(qlist)) {
    fail("likes file and Q-list file have different lengths", 5)
}

best <- which.max(ll)
Qb <- read_q_matrix(qlist[best])

results <- vector("list", length(ll))
for (i in seq_along(ll)) {
    Qt <- read_q_matrix(qlist[i])
    Qt <- get_fast(Qt, Qb)
    diff <- abs(Qt - Qb)
    diffsum <- colSums(diff)
    results[[i]] <- data.frame(
        run = i,
        ll = ll[i],
        lldiff = ll[i] - ll[best],
        mean = mean(diffsum),
        maxdiffsum = max(diffsum),
        rmse = sqrt(mean(colSums((Qt - Qb)^2))),
        q95 = as.numeric(quantile(diffsum, 0.95)),
        q975 = as.numeric(quantile(diffsum, 0.975)),
        max = max(diff),
        conv = max(diff) < threshold
    )
}

res <- do.call(rbind, results)
conv_count <- sum(res$conv)

if (count_only) {
    cat(conv_count, "\n", sep = "")
} else {
    write.table(res, file = stdout(), sep = "\t", quote = FALSE, row.names = FALSE)
}
