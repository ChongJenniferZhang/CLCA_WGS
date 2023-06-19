# Please run this script from the top directory
if (basename(getwd()) != "CLCA_WGS") {
  stop("Please run from top level directory, CLCA_WGS")
}

# This script is designed to move only important SigProfierExtractor output
# files to upper level directory. The rest of the raw output will be git ignored
# as their paths are too long and may cause error when using git for version control

sigpro_dir_dbs78 <-
  "signature/extraction/raw_output/DBS78/SigProfilerExtractor"

paths_to_copy <-
  c("DBS78/DBS78_selection_plot.pdf",
    "DBS78/Suggested_Solution/DBS78_De-Novo_Solution/Signatures/DBS78_De-Novo_Signatures.txt",
    "DBS78/Suggested_Solution/DBS78_De-Novo_Solution/Signatures/DBS_78_plots_DBS78_De-Novo.pdf",
    "DBS78/Suggested_Solution/COSMIC_DBS78_Decomposed_Solution/DBS78_Decomposition_Plots.pdf",
    "DBS78/Suggested_Solution/DBS78_De-Novo_Solution/Activities/DBS78_De-Novo_Activities_refit.txt")

sapply(X = paths_to_copy, FUN = function(path_to_copy) {
  file_to_move <- file.path(sigpro_dir_dbs78, path_to_copy)
  file.copy(
    from = file_to_move, to = sigpro_dir_dbs78,
    overwrite = TRUE, copy.date = TRUE
  )
})