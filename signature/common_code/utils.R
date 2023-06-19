
classify_sample_sig <- function(mut_type) {
  load("signature/attribution/raw_output/SBS96/sbs_h2_h8_samples.Rdata")
  catalog_file <- file.path(
    "signature/attribution/input", mut_type,
    paste0("clca_catalog_", mut_type, ".csv")
  )
  sig_file <- file.path(
    "signature/attribution/input", mut_type,
    paste0("clca_signature_", mut_type, ".csv")
  )

  catalog <- ICAMS::ReadCatalog(file = catalog_file)
  sig <- ICAMS::ReadCatalog(
    file = sig_file,
    catalog.type = "counts.signature"
  )

  catalogs <- sigs <- list()

  if (mut_type == "SBS96") {
    sbs_h2_h8_both_catalog <- catalog[, sbs_h2_h8_both]
    sbs_h8_only_catalog <- catalog[, sbs_h8_only]
    sbs_h2_only_catalog <- catalog[, sbs_h2_only]
    sbs_h2_h8_neither <-
      setdiff(colnames(catalog), c(sbs_h2_h8_both, sbs_h2_only, sbs_h8_only))
    sbs_h2_h8_neither_catalog <- catalog[, sbs_h2_h8_neither]

    sbs_h2_h8_both_sig <- sig
    sbs_h8_only_sig <- sig[, colnames(sig) != "SBS_H2"]
    sbs_h2_only_sig <- sig[, colnames(sig) != "SBS_H8"]
    sbs_h2_h8_neither_sig <- sig[, !colnames(sig) %in% c("SBS_H2", "SBS_H8")]

    catalogs[["sbs_h2_h8_both"]] <- sbs_h2_h8_both_catalog
    catalogs[["sbs_h8_only"]] <- sbs_h8_only_catalog
    catalogs[["sbs_h2_only"]] <- sbs_h2_only_catalog
    catalogs[["sbs_h2_h8_neither"]] <- sbs_h2_h8_neither_catalog

    sigs[["sbs_h2_h8_both"]] <- sbs_h2_h8_both_sig
    sigs[["sbs_h8_only"]] <- sbs_h8_only_sig
    sigs[["sbs_h2_only"]] <- sbs_h2_only_sig
    sigs[["sbs_h2_h8_neither"]] <- sbs_h2_h8_neither_sig
  } else {
    # Use the SBS96 attribution result to classify samples and signatures
    sbs96_exposure_file <- "signature/attribution/summary/clca_exposure_SBS96.csv"
    if (!file.exists(sbs96_exposure_file)) {
      stop("Please do signature attribution for SBS96 first and save the exposure file")
    }

    sbs96_exposure <- mSigTools::read_exposure(file = sbs96_exposure_file)
    sbs_h2_exposure <- sbs96_exposure["SBS_H2", ]
    sbs_h2_sample <- names(sbs_h2_exposure[sbs_h2_exposure > 0])
    sbs_h3_exposure <- sbs96_exposure["SBS_H3", ]
    sbs_h3_sample <- names(sbs_h3_exposure[sbs_h3_exposure > 0])

    if (mut_type == "DBS78") {
      # Both SBS_H2 and DBS_H2 are aristolochic acid (AA) signatures
      dbs_h2_catalog <- catalog[, sbs_h2_sample]
      non_sbs_h2_sample <- setdiff(colnames(catalog), sbs_h2_sample)
      non_dbs_h2_catalog <- catalog[, non_sbs_h2_sample]

      dbs_h2_sig <- sig
      non_dbs_h2_sig <- sig[, colnames(sig) != "DBS_H2"]

      catalogs[["dbs_h2"]] <- dbs_h2_catalog
      catalogs[["non_dbs_h2"]] <- non_dbs_h2_catalog

      sigs[["dbs_h2"]] <- dbs_h2_sig
      sigs[["non_dbs_h2"]] <- non_dbs_h2_sig
    } else if (mut_type == "ID83") {
      sbs_h2_h3_both_sample <- intersect(sbs_h2_sample, sbs_h3_sample)
      sbs_h2_only_sample <- setdiff(sbs_h2_sample, sbs_h2_h3_both_sample)
      sbs_h3_only_sample <- setdiff(sbs_h3_sample, sbs_h2_h3_both_sample)
      sbs_h2_h3_neither_sample <-
        setdiff(
          colnames(catalog),
          c(sbs_h2_h3_both_sample, sbs_h2_only_sample, sbs_h3_only_sample)
        )

      # Both SBS_H2 and ID_H3 are aristolochic acid (AA) signatures
      # Both SBS_H3 and ID_H8 are aflatoxin signatures
      id_h3_h8_both_catalog <- catalog[, sbs_h2_h3_both_sample]
      id_h3_only_catalog <- catalog[, sbs_h2_only_sample]
      id_h8_only_catalog <- catalog[, sbs_h3_only_sample]
      id_h3_h8_neither_catalog <- catalog[, sbs_h2_h3_neither_sample]

      id_h3_h8_both_sig <- sig
      id_h3_only_sig <- sig[, colnames(sig) != "ID_H8"]
      id_h8_only_sig <- sig[, colnames(sig) != "ID_H3"]
      id_h3_h8_neither_sig <- sig[, !colnames(sig) %in% c("ID_H3", "ID_H8")]

      catalogs[["id_h3_h8_both"]] <- id_h3_h8_both_catalog
      catalogs[["id_h3_only"]] <- id_h3_only_catalog
      catalogs[["id_h8_only"]] <- id_h8_only_catalog
      catalogs[["id_h3_h8_neither"]] <- id_h3_h8_neither_catalog

      sigs[["id_h3_h8_both"]] <- id_h3_h8_both_sig
      sigs[["id_h3_only"]] <- id_h3_only_sig
      sigs[["id_h8_only"]] <- id_h8_only_sig
      sigs[["id_h3_h8_neither"]] <- id_h3_h8_neither_sig
    }
  }

  return(list(catalogs = catalogs, sigs = sigs))
}

check_sort_exposure <- function(exposure, catalog) {
  names1 <- colnames(exposure)
  names2 <- colnames(catalog)
  if (setequal(names1, names2)) {
    return(exposure[, names2, drop = FALSE])
  } else {
    zero_mutation_samples <- setdiff(names2, names1)
    zero_exposure_matrix <-
      matrix(
        data = 0, nrow = nrow(exposure), ncol = length(zero_mutation_samples),
        dimnames = list(rownames(exposure), zero_mutation_samples)
      )
    tmp <- cbind(exposure, zero_exposure_matrix)
    return(tmp[, names2, drop = FALSE])
  }
}

combine_write_exposure <- function(mut_type) {
  output_dir <-
    file.path("signature/attribution/raw_output", mut_type, "mSigAct")
  rds_files <-
    list.files(
      path = output_dir, pattern = "\\.Rds$", full.names = TRUE,
      recursive = TRUE
    )

  exposure_list <- lapply(rds_files, FUN = function(rds_file) {
    retval <- readRDS(rds_file)
    return(retval$proposed.assignment)
  })

  exposure <- mSigAct:::MergeListOfExposures(exposure_list)

  input_file <-
    file.path(
      "signature/attribution/input", mut_type,
      paste0("clca_catalog_", mut_type, ".csv")
    )

  catalog <- ICAMS::ReadCatalog(input_file)

  # Sort exposure according to the original order of sample names
  exposure2 <- check_sort_exposure(exposure = exposure, catalog = catalog)

  output_dir <- "signature/attribution/summary"
  if (!dir.exists(output_dir)) {
    dir.create(path = output_dir, recursive = TRUE)
  }

  output_file <- paste0("clca_exposure_", mut_type, ".csv")
  mSigTools::write_exposure(
    exposure = exposure2,
    file = file.path(output_dir, output_file)
  )
}
