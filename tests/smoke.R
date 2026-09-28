#!/usr/bin/env Rscript
# End-to-end smoke test: write a small low-rank gene x cell matrix as a TENx h5
# (what FEAT emits), run pca.R once per solver, check the output contract and that
# irlba/random recover exact's leading PCs. Run from the module root:
#   Rscript tests/smoke.R
suppressPackageStartupMessages({
  library(Matrix)
  library(HDF5Array)
  library(data.table)
})

set.seed(1)
n_genes <- 300; n_cells <- 200; k <- 10
signal <- matrix(rnorm(n_genes * 5), n_genes) %*% matrix(rnorm(5 * n_cells), 5) * 3
X <- pmax(signal + matrix(rnorm(n_genes * n_cells), n_genes), 0)  # non-negative, ~half zeros
X <- as(X, "dgCMatrix")
dimnames(X) <- list(sprintf("gene%03d", seq_len(n_genes)), sprintf("cell%03d", seq_len(n_cells)))

tmp <- tempfile("pcabs"); dir.create(tmp)
h5 <- file.path(tmp, "t_normalized_selected.h5")
invisible(writeTENxMatrix(X, h5, group = "matrix"))

run <- function(solver) {
  out <- file.path(tmp, solver)
  status <- system2("Rscript", c("pca.R", "--output_dir", out, "--name", "t",
                                 "--normalized_selected_h5", h5, "--solver", solver,
                                 "--n_components", k, "--random_seed", 42,
                                 if (solver == "random") c("--power_iters", 7)),
                    stdout = FALSE, stderr = FALSE)
  stopifnot(status == 0)
  emb <- fread(file.path(out, "t_embedding.tsv"))
  load <- fread(file.path(out, "t_loadings.tsv"))
  # contract: one row per cell/gene in input order, PC1..PCk, finite
  stopifnot(identical(names(emb), c("cell_id", paste0("PC", 1:k))),
            identical(emb$cell_id, colnames(X)),
            identical(names(load), c("gene", paste0("PC", 1:k))),
            identical(load$gene, rownames(X)),
            all(is.finite(as.matrix(emb[, -1]))), all(is.finite(as.matrix(load[, -1]))))
  as.matrix(emb[, -1])
}

ref <- run("exact")
for (solver in c("irlba", "random")) {
  e <- run(solver)
  # the 5 signal PCs must agree with exact up to sign
  agree <- abs(diag(cor(ref[, 1:5], e[, 1:5])))
  stopifnot(all(agree > 0.99))
  cat(sprintf("%-6s min |cor| vs exact over PC1-5: %.4f\n", solver, min(agree)))
}
unlink(tmp, recursive = TRUE)
cat("OK\n")
