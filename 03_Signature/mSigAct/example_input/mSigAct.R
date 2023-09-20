.libPaths(c("/usr/local/lib/R/site-library", # R library path inside container
            .libPaths()))

# args have the arguments passed to Rscript from bash script
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Please provide the random seed as an argument")
}

seed <- as.numeric(args[1])

library(ICAMS)
library(mSigAct)
library(cosmicsig)

input_file <- "example_input/catalog.csv"
input_catalog <- ICAMS::ReadCatalog(file = input_file)
test_catalog <- input_catalog[, 1:2, drop = FALSE]
sigs <- cosmicsig::COSMIC_v3.2$signature$GRCh37$ID
sig_prop <- mSigAct::ExposureProportions(mutation.type = "ID", 
                                         cancer.type = "Biliary-AdenoCA")
sigs_to_use <- sigs[, names(sig_prop), drop = FALSE]

output_home <- "example_output/mSigAct"

retval <-
  mSigAct:: MAPAssignActivity(spectra = test_catalog, 
                              sigs = sigs_to_use, 
                              sigs.presence.prop = sig_prop,
                              output.dir = output_home, 
                              max.level = ncol(sigs_to_use) - 1,
                              p.thresh = 0.05 / ncol(sigs_to_use), 
                              num.parallel.samples = 2, 
                              mc.cores.per.sample = 5, 
                              seed = seed)
