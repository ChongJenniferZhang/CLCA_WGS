# This script is designed to run in a singularity/apptainer container,
# it is suggested to use a bash script to call this R script

# We need to tell R to look for packages first inside the container
# but also to be able to look outside the container (the "host"
# system). This works only if the host OS is the same as
# the container OS, or if the libraries in the host system do not
# have compiled code. In the last case, you would have to
# convert the container to a sandbox, add the necessary libraries,
# and then convert back to a .sif.
.libPaths(c(
  "/usr/local/lib/R/site-library", # R library path inside the container
  .libPaths()
))

# args have the arguments passed to Rscript from bash script
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Please provide the project directory as an argument")
}
project_dir <- args[1]

library(mSigAct)
library(ICAMS)
library(mSigTools)

setwd(project_dir)

source("signature/common_code/utils.R")

mut_type <- "DBS78"
combine_write_exposure(mut_type)
