library(CHASM)

# Example package calls derived from the manual scripts in
# scripts/v7_Wavelet_transform/manual_tests/.

chromosomes <- paste0("chr", c(1:22, "X", "Y"))

# Wavelet pipeline:
# chrom_depth_per_cell <- normalize_depth(df.depth = cn_bin.wide, bins = positions)
# wt <- wavelet_transform(chrom_depth_per_cell, bins = positions, chromosomes = chromosomes)
# rpca <- robust_pca(wt$mat.wavelet.transform)
# sparse_signal <- inv_wavelet_transform(rpca$Sparse_Signal, wt$chrom.informed.wavelet)
# segment_output <- segment_residuals(sparse_signal)
# cn_states <- assign_cn_state(chrom_depth_per_cell, sparse_signal, segment_output)

# Chromosome-level negative-binomial check:
# cn_state_nb <- assign_cn_state.chrom(read_depth.chrom, positions)
