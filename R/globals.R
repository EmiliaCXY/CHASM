utils::globalVariables(c(
  "ID",
  "CN",
  "CN_state.binom",
  "CN_state_raw",
  "barcode",
  "beta_i_c",
  "bin",
  "calledCNA",
  "celltype",
  "chrom",
  "chrom_arm.adjust",
  "chrom_name",
  "depth_bin",
  "depth_bin_new",
  "expected_count",
  "identifier",
  "lib_size",
  "n_segments",
  "observed_count",
  "p_val_adjust",
  "read_depth",
  "segment_output_ID",
  "start",
  "total_in_group",
  "weight"
))

# data.table non-standard evaluation in R/second_round.R
utils::globalVariables(c("chrom", "variable", "segment_id", "cn", "ID", "nseg", "dir", "n", "i.n", "hit", "run"))
.datatable.aware <- TRUE
