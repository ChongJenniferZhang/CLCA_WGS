# This script is designed to run in a singularity/apptainer container,
# it is suggested to use a bash script to call this R script

# We need to tell R to look for packages first inside the container
# but also to be able to look outside the container (the "host"
# system). This works only if the host OS is the same as
# the container OS, or if the libraries in the host system do not
# have compiled code. In the last case, you would have to
# convert the container to a sandbox, add the necessary libraries,
# and then convert back to a .sif.
.libPaths(c("/usr/local/lib/R/site-library", # R library path inside the container
            .libPaths()))

# args have the arguments passed to Rscript from bash script
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Please provide the project directory as an argument")
}
project_dir <- args[1]

library(ICAMS)
library(mSigAct)

setwd(project_dir)

source("signature/common_code/utils.R")

mut_type <- "ID83"

inputs <- classify_sample_sig(mut_type)
output_home <- file.path("signature/attribution/raw_output", mut_type)

message("Start running mSigAct")
message("mSigAct version: ")
print(packageVersion("mSigAct"))

time_used <- system.time(expr = {

for (type in names(inputs$catalogs)) {
  catalog <- inputs$catalogs[[type]]
  sig <- inputs$sigs[[type]]
  output_dir <- file.path(output_home, "raw", type)
  
  retval <-
    mSigAct::SparseAssignActivity(spectra = catalog, 
                                  sigs = sig, 
                                  output.dir = output_dir, 
                                  max.level = ncol(sig) - 1, 
                                  p.thresh = 0.05 / ncol(sig), 
                                  num.parallel.samples = 12, 
                                  mc.cores.per.sample = 8, 
                                  seed = 5189, 
                                  max.subsets = .Machine$double.xmax, 
                                  drop.low.mut.samples = FALSE)
  output_folder <- file.path(output_home, "mSigAct")
  if (!dir.exists(output_folder)) {
    dir.create(path = output_folder, recursive = TRUE)
  }
  saveRDS(retval, 
          file = file.path(output_folder, paste0("retval_", type, ".Rds")))
}
})

saveRDS(time_used, file = file.path(output_home, "time_used.Rds"))

message("Time used: ")
print(time_used)
message("Finished running mSigAct")
