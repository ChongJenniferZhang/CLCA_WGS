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

library(ICAMS)
library(mSigAct)

input_catalog <-
  file.path(project_dir, "signature/attribution/input/SBS96/clca_catalog_SBS96.csv")
input_sig <-
  file.path(project_dir, "signature/attribution/input/SBS96/clca_signature_SBS96.csv")
sbs_catalog <- ICAMS::ReadCatalog(file = input_catalog)
sbs_sig <- ICAMS::ReadCatalog(
  file = input_sig,
  catalog.type = "counts.signature"
)

my_opts <- mSigAct::DefaultManyOpts(likelihood.dist = "neg.binom")
my_opts$nbinom.size <- 4

# Signature presence test for SBS_H2 (SBS22)
sbs_h2_retval <-
  mSigAct::SignaturePresenceTest(
    spectra = sbs_catalog,
    sigs = sbs_sig,
    target.sig.index = "SBS_H2",
    m.opts = my_opts,
    seed = 5189,
    mc.cores = 100
  )
sbs_h2_p_values <- sapply(sbs_h2_retval, FUN = "[[", 4)
sbs_h2_q_values <- stats::p.adjust(p = sbs_h2_p_values, method = "BH")
sbs_h2_samples <- names(sbs_h2_q_values[sbs_h2_q_values < 0.05])

# Signature presence test for SBS_H8
sbs_h8_retval <-
  mSigAct::SignaturePresenceTest(
    spectra = sbs_catalog,
    sigs = sbs_sig,
    target.sig.index = "SBS_H8",
    m.opts = my_opts,
    seed = 5189,
    mc.cores = 100
  )
sbs_h8_p_values <- sapply(sbs_h8_retval, FUN = "[[", 4)
sbs_h8_q_values <- stats::p.adjust(p = sbs_h8_p_values, method = "BH")
sbs_h8_samples <- names(sbs_h8_q_values[sbs_h8_q_values < 0.05])

# Find out samples which have both SBS_H2 and SBS_H8
sbs_h2_h8_both <- intersect(sbs_h2_samples, sbs_h8_samples)

sbs_h2_only <- setdiff(sbs_h2_samples, sbs_h2_h8_both)
sbs_h8_only <- setdiff(sbs_h8_samples, sbs_h2_h8_both)

output_dir <- file.path(project_dir, "signature/attribution/raw_output/SBS96")
if (!dir.exists(output_dir)) {
  dir.create(path = output_dir, recursive = TRUE)
}
save(sbs_h2_h8_both, sbs_h2_only, sbs_h8_only,
  file = file.path(output_dir, "sbs_h2_h8_samples.Rdata")
)
