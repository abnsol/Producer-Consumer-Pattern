#!/usr/bin/env Rscript
# Installs required R packages for the pipeline.

options(warn = 1)

# CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

install_if_missing <- function(pkgs) {
  to_install <- pkgs[!(pkgs %in% rownames(installed.packages()))]
  if (length(to_install)) {
    install.packages(to_install, dependencies = TRUE)
  }
}

# Core build deps often needed by packages
install_if_missing(c("Rcpp", "RcppArmadillo", "Matrix", "data.table"))

# susieR is on CRAN
install_if_missing(c("susieR"))

# MungeSumstats is on Bioconductor
if (!"BiocManager" %in% rownames(installed.packages())) {
  install.packages("BiocManager")
}
BiocManager::install("MungeSumstats", update = FALSE, ask = FALSE)

# Optional helpful packages
install_if_missing(c("devtools", "remotes"))