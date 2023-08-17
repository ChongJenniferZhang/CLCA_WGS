# Please run this script from the top directory
if (basename(getwd()) != "CLCA_WGS") {
  stop("Please run from top level directory, CLCA_WGS")
}

# This script is designed to move only important SigProfierExtractor output
# files to upper level directory. The rest of the raw output will be git ignored
# as their paths are too long and may cause error when using git for version control

sigpro_dir_id83 <-
  "signature/extraction/raw_output/ID83/SigProfilerExtractor"

paths_to_copy <-
  c("ID83/ID83_selection_plot.pdf",
    "ID83/Suggested_Solution/ID83_De-Novo_Solution/Signatures/ID83_De-Novo_Signatures.txt",
    "ID83/Suggested_Solution/ID83_De-Novo_Solution/Signatures/ID_83_plots_ID83_De-Novo.pdf",
    "ID83/Suggested_Solution/COSMIC_ID83_Decomposed_Solution/ID83_Decomposition_Plots.pdf",
    "ID83/Suggested_Solution/ID83_De-Novo_Solution/Activities/ID83_De-Novo_Activities_refit.txt")

sapply(X = paths_to_copy, FUN = function(path_to_copy) {
  file_to_move <- file.path(sigpro_dir_id83, path_to_copy)
  file.copy(
    from = file_to_move, to = sigpro_dir_id83,
    overwrite = TRUE, copy.date = TRUE
  )
})