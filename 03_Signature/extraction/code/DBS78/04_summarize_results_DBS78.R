# Please run this script from the top directory
if (basename(getwd()) != "CLCA_WGS") {
  stop("Please run from top level directory, CLCA_WGS")
}

library(mSigTools) # remotes::install_github(repo = "Rozen-Lab/mSigTools", ref = "v1.0.6-branch")
library(ICAMS) # remotes::install_github(repo = "steverozen/ICAMS", ref = "v3.0.6-branch")
library(gtools)

msighdp_sig_file <-
  "signature/extraction/raw_output/DBS78/mSigHdp/extracted.signatures.csv"
msighdp_sig <-
  ICAMS::ReadCatalog(file = msighdp_sig_file, catalog.type = "counts.signature")

# Remove potential mSigHdp-extracted signature
indices <- grep(pattern = "potential", x = colnames(msighdp_sig))
msighdp_sig2 <- msighdp_sig[, -indices, drop = FALSE]

sigpro_sig_file <-
  "signature/extraction/raw_output/DBS78/SigProfilerExtractor/DBS78_De-Novo_Signatures.txt"
sigpro_sig <-
  ICAMS::ReadCatalog(file = sigpro_sig_file, catalog.type = "counts.signature")

# Check how many mSigHdp extracted signatures can be matched to SigProfiler
# extracted signatures
match_retval <- mSigTools::TP_FP_FN_avg_sim(
  extracted.sigs = sigpro_sig,
  reference.sigs = msighdp_sig2,
  similarity.cutoff = -1
)
match_table <- match_retval$table

# Find mSigHdp-extracted signature that has a cosine similarity ≥ 0.90 with a
# SigProfiler-extracted signature
ids_to_include <- match_table[match_table$sim >= 0.9, ]$ref.sig

ids_to_check <- setdiff(colnames(msighdp_sig2), ids_to_include)
match_table_to_check <-
  match_table[match_table$ref.sig %in% ids_to_check, ]

# Check whether the mSigHdp-extracted signature can be reconstructed by multiple
# SigProfiler-extracted signatures (reconstruction cosine similarity≥0.90)

# Column sums of SigProfiler-extracted signatures are not exactly 1 (though
# difference is very small). We need to make sure the column sums are exactly 1
# in order to use function mSigTools::best_reconstruction_QP

all.equal.numeric(colSums(sigpro_sig),
  rep(1, ncol(sigpro_sig)),
  check.attributes = FALSE
)

# Normalize SigProfiler-extracted signatures to make sure column sums are exactly 1
sigpro_sig2 <-
  sweep(x = sigpro_sig, MARGIN = 2, STATS = colSums(sigpro_sig), FUN = "/")

retval <- lapply(X = ids_to_check, FUN = function(id_to_check) {
  msighdp_sig_to_check <- msighdp_sig[, id_to_check, drop = FALSE]
  reconstruction <-
    mSigTools::best_reconstruction_QP(
      target.sig = msighdp_sig_to_check,
      sig.universe = sigpro_sig2,
      max.subset.size = 3
    )
  return(reconstruction)
})
names(retval) <- ids_to_check
cosines <- sapply(X = retval, FUN = "[", 2)
all(cosines > 0.9)

ids_final <- c(ids_to_include, ids_to_check)
msighdp_sig_final <- msighdp_sig[, colnames(msighdp_sig) %in% ids_final]

# Rename and reorder the msighdp extracted signatures
colnames(msighdp_sig_final) <- c("DBS_H2", "DBS_H3", "DBS_H1")

msighdp_sig_final <-
  msighdp_sig_final[, gtools::mixedsort(colnames(msighdp_sig_final))]

ICAMS::WriteCatalog(
  catalog = msighdp_sig_final,
  file = "signature/extraction/summary/clca_signature_DBS78.csv"
)
