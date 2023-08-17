# Please run this script from the top directory
if (basename(getwd()) != "CLCA_WGS") {
  stop("Please run from top level directory, CLCA_WGS")
}

library(mSigTools) # remotes::install_github(repo = "Rozen-Lab/mSigTools", ref = "v1.0.6-branch")
library(ICAMS) # remotes::install_github(repo = "steverozen/ICAMS", ref = "v3.0.6-branch")
library(gtools)

msighdp_sig_file <-
  "signature/extraction/raw_output/SBS96/mSigHdp/extracted.signatures.csv"
msighdp_sig <-
  ICAMS::ReadCatalog(file = msighdp_sig_file, catalog.type = "counts.signature")

sigpro_sig_file <-
  "signature/extraction/raw_output/SBS96/SigProfilerExtractor/SBS96_De-Novo_Signatures.txt"
sigpro_sig <-
  ICAMS::ReadCatalog(file = sigpro_sig_file, catalog.type = "counts.signature")

# Check how many mSigHdp extracted signatures can be matched to SigProfiler
# extracted signatures
match_retval <- mSigTools::TP_FP_FN_avg_sim(
  extracted.sigs = sigpro_sig,
  reference.sigs = msighdp_sig,
  similarity.cutoff = -1
)
match_table <- match_retval$table

# Discard those mSigHdp extracted signatures that cannot be matched to any
# SigProfiler extracted signature
ids_to_discard <- setdiff(colnames(msighdp_sig), match_table$ref.sig)
msighdp_sig2 <- msighdp_sig[, !colnames(msighdp_sig) %in% ids_to_discard]

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

ids_union <- c(ids_to_include, ids_to_check)
msighdp_sig_union <- msighdp_sig[, ids_union]
msighdp_sig_union_file <-
  "signature/extraction/raw_output/SBS96/msighdp_sig_union.pdf"

ICAMS::PlotCatalogToPdf(
  catalog = msighdp_sig_union,
  file = msighdp_sig_union_file
)

# By manual inspection, we thought hdp.17 is a split signature of hdp.7
# with T>C and T>G peaks missing. So we do not include hdp.17 in the final
# signature set
ids_final <- ids_union[ids_union != "hdp.17"]
unlink(x = msighdp_sig_union_file, recursive = TRUE)

msighdp_sig_final <- msighdp_sig[, colnames(msighdp_sig) %in% ids_final]

# Rename and reorder the msighdp extracted signatures
colnames(msighdp_sig_final) <-
  c(
    "SBS_H2", "SBS_H3", "SBS_H5", "SBS_H4", "SBS_H6", "SBS_H7",
    "SBS_H8", "SBS_H9", "SBS_H10", "SBS_H11", "SBS_H12", "SBS_H13",
    "SBS_H14", "SBS_H16", "SBS_H15", "SBS_H1", "SBS_H17"
  )

msighdp_sig_final <-
  msighdp_sig_final[, gtools::mixedsort(colnames(msighdp_sig_final))]

ICAMS::WriteCatalog(
  catalog = msighdp_sig_final,
  file = "signature/extraction/summary/clca_signature_SBS96.csv"
)

