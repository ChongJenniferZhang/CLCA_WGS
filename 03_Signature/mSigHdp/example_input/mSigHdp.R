.libPaths(c("/usr/local/lib/R/site-library", 
            .libPaths()))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Please provide the random seed as an argument")
}

seed <- as.numeric(args[1])

library(ICAMS)
library(mSigHdp)
library(cosmicsig)

input_file <- "example_input/catalog.csv"
input_catalog <- ICAMS::ReadCatalog(file = input_file)

output_home <- "example_output"

if (!dir.exists(output_home)) {
  if (!dir.create(path = output_home, recursive = TRUE)) {
    stop("unable to create ", output_home)
  }
}

setwd(output_home)

K_guess <- 15

mSigHdp::RunHdpxParallel(input.catalog      = input_catalog, 
                         ground.truth.sig   = cosmicsig::COSMIC_v3.2$signature$GRCh37$ID, 
                         out.dir            = output_home, 
                         CPU.cores          = 2, 
                         num.child.process  = 2, 
                         seedNumber         = seed, 
                         K.guess            = K_guess,
                         burnin.checkpoint  = TRUE, 
                         burnin             = 100,  
                         burnin.multiplier  = 2,
                         post.n             = 5, 
                         post.space         = 5,
                         multi.types        = FALSE, 
                         overwrite          = TRUE, 
                         gamma.alpha        = 1,
                         gamma.beta         = 20,
                         cos.merge          = 0.90,
                         confident.prop     = 0.6) 