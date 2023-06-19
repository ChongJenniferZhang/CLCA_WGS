# Please run this script from the top directory
if (basename(getwd()) != "CLCA_WGS") {
  stop("Please run from top level directory, CLCA_WGS")
}

# This script is designed to move only important SigProfierExtractor output
# files to upper level directory. The rest of the raw output will be git ignored
# as their paths are too long and may cause error when using git for version control

sigpro_dir_sbs96 <-
  "signature/extraction/raw_output/SBS96/SigProfilerExtractor"

paths_to_copy <-
  c("SBS96/SBS96_selection_plot.pdf",
    "SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Signatures/SBS96_De-Novo_Signatures.txt",
    "SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Signatures/SBS_96_plots_SBS96_De-Novo.pdf",
    "SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/SBS96_Decomposition_Plots.pdf",
    "SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Activities/SBS96_De-Novo_Activities_refit.txt")

sapply(X = paths_to_copy, FUN = function(path_to_copy) {
  file_to_move <- file.path(sigpro_dir_sbs96, path_to_copy)
  file.copy(
    from = file_to_move, to = sigpro_dir_sbs96,
    overwrite = TRUE, copy.date = TRUE
  )
})
