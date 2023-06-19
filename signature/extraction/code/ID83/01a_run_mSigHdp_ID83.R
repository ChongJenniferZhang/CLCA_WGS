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
library(mSigHdp)
library(cosmicsig)

input_file <-
  file.path(project_dir, "signature/extraction/input/ID83/clca_catalog_ID83.csv")
indel_catalog <- ICAMS::ReadCatalog(file = input_file)

output_dir <- file.path(project_dir, "signature/extraction/raw_output/ID83/mSigHdp/")

if (!dir.exists(output_dir)) {
  dir.create(path = output_dir, recursive = TRUE)
}

checkpoint_dir <- file.path(output_dir, "checkpoint")
if (!dir.exists(checkpoint_dir)) {
  dir.create(path = checkpoint_dir, recursive = TRUE)
}
# Set working directory to checkpoint_dir directory so that the checkpoint Rdata
# files will be saved there
setwd(checkpoint_dir)

message("Start running mSigHdp")
message("Using seed 1234")
message("Using K.guess 10")

message("mSigHdp version: ")
print(packageVersion("mSigHdp"))

num_jobs <- 20 #Defines how many independent posterior chains were initiated. 20 is recommended.

time_used <- system.time(expr = {
  mSigHdp::RunHdpxParallel(input.catalog = indel_catalog, 
                           ground.truth.sig  = cosmicsig::COSMIC_v3.2$signature$GRCh37$ID, #use COSMIC signatures as ground truth signatures
                           out.dir            = output_dir, # the directory of folder to store the results
                           CPU.cores          = num_jobs, #number of cores for parallel running
                           num.child.process  = num_jobs, #number of chains for parallel running
                           seedNumber         = 1234, # for randomization
                           K.guess            = 10, #estimated initiating number of signatures
                           burnin.checkpoint  = TRUE, #checkpoint for burnin iterations
                           burnin        = 5000,  
                           burnin.multiplier  = 3,
                           post.n             = 200, #number of posterior samples to be collected
                           post.space         = 100, #number of iterations between posterior samples
                           multi.types        = FALSE, # build a two-layer HDP structure
                           overwrite          = TRUE, #overwrite if the folder has old results
                           gamma.alpha = 1,#gamma distribution parameter (shape) for concentration parameter
                           gamma.beta  = 50,#gamma distribution parameter (rate) for concentration parameter
                           cos.merge = 0.88,#a cutoff to merge similar signatures (cosine similarity >= 0.9)
                           confident.prop = 0.9) #a signature is reported if it is found in no less than 60% of posterior samples
}
)

saveRDS(time_used, file = file.path(output_dir, "time_used.Rds"))
message("Time used: ")
print(time_used)
message("Finished running mSigHdp")
