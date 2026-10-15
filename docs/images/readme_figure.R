# Makes docs/images/Demo2_admixture_K2-4.png (README figure). Run from the repository root after:
#   for K in 2 3 4; do admixer -j8 --seed=1 --conv 3 -o Demo/Results/Admixer/Demo2fig Demo/Data/Demo2input.gz $K; done
#   Rscript docs/images/readme_figure.R Demo/Results/Admixer/Demo2fig docs/images/Demo2_admixture_K2-4.png
args <- commandArgs(TRUE); prefix <- args[1]; out <- args[2]
pop <- read.table("Demo/Data/Demo2pop.info", as.is = TRUE)[, 1]
popOrder <- c("YRI", "ASW", "CEU", "MXL", "CHB")
Ks <- 2:4
cols <- c("#D55E00", "#0072B2", "#009E73", "#E69F00", "#CC79A7")
Q <- lapply(Ks, function(K) as.matrix(read.table(paste0(prefix, ".", K, ".Q"))))
# match columns of each K to the previous K so colours stay the same
for (i in 2:length(Q)) {
  prev <- Q[[i - 1]]; cur <- Q[[i]]; K <- ncol(cur); left <- 1:K; ordc <- integer(0)
  for (j in 1:ncol(prev)) {
    best <- left[which.max(sapply(left, function(l) sum(abs(prev[, j] - cur[, l]) < 0.1)))]
    ordc <- c(ordc, best); left <- setdiff(left, best)
  }
  Q[[i]] <- cur[, c(ordc, left)]
}
last <- Q[[3]]
ord <- unlist(lapply(popOrder, function(p) { i <- which(pop == p); i[order(-last[i, which.max(colMeans(last[i, , drop = FALSE]))])] }))
png(out, width = 2000, height = 800, res = 200)
par(mfrow = c(length(Ks), 1), mar = c(0.4, 4, 0.4, 1), oma = c(2.2, 0, 2, 0))
for (i in seq_along(Ks)) {
  barplot(t(Q[[i]][ord, ]), col = cols, space = 0, border = NA, axes = FALSE, xaxs = "i")
  mtext(paste0("K = ", Ks[i]), side = 2, line = 1, las = 1, cex = 0.8, adj = 0.75)
  b <- cumsum(table(factor(pop[ord], levels = popOrder)))
  abline(v = b[-length(b)], col = "white", lwd = 2)
}
mid <- b - table(factor(pop[ord], levels = popOrder)) / 2
mtext(popOrder, side = 1, at = mid, line = 0.8, cex = 0.8)
mtext("Admixture proportions estimated from genotype likelihoods (NGSadmix model, admixer)", side = 3, outer = TRUE, line = 0.5, cex = 0.85)
invisible(dev.off())
