# CHASM

`CHASM` packages the wavelet-transform, segmentation, and copy-number inference
functions from the `scATACcnv` project into reusable R functions.

## Installation

Install the required R packages:

```r
install.packages(c(
  "DNAcopy",
  "dplyr",
  "GenomicRanges",
  "magrittr",
  "stringr",
  "tidyr",
  "wavethresh",
  "rsvd"
))
```

Then install `CHASM` from the package directory:

```r
install.packages("/path/to/scATACcnv/package", repos = NULL, type = "source")
```

For development use, you can load the package without installing it:

```r
pkgload::load_all("/path/to/scATACcnv/package")
```

## Example Data

The package bundles an example read-depth matrix at
`inst/extdata/Read_depth_matrix_chr1_full_2_5pct_Rep1.rds`. Load it with:

```r
library(CHASM)

example_rds <- example_read_depth_path()
read_depth <- readRDS(example_rds)
read_depth$barcode <- rownames(read_depth)
read_depth <- read_depth[, c("barcode", setdiff(colnames(read_depth), "barcode"))]
```

This produces a wavelet-workflow input table with one row per cell, a
`barcode` column, and one column per genomic bin.

To discover other bundled examples:

```r
list_example_data()

read_depth_chrom <- readRDS(
  example_data_path("Read_depth_matrix_chr1_full_2_5pct_Rep1_chrom.rds")
)
```

## Input Data

`CHASM` currently supports two related workflows.

### 1. Wavelet-based segment-level CN calling

The wavelet workflow expects a read-depth matrix with:

- one row per cell
- row names corresponding to cell barcodes
- one `barcode` column
- one column per genomic bin

Bin names should follow the package convention, for example:

```r
"chr1_1_2000000_p"
"chr1_2000001_4000000_q"
```

### 2. Chromosome-level negative-binomial CN calling

The chromosome-level workflow expects a read-depth table with:

- one row per cell
- one `ID` column
- one `celltype` column
- one column per chromosome-level feature listed in `positions`

## Usage

### Wavelet workflow

```r
library(CHASM)

chromosomes <- paste0("chr", c(1:22, "X", "Y"))

bins <- setdiff(colnames(read_depth), "barcode")
chrom_depth <- normalize_depth(read_depth, bins, chromosomes)
wt <- wavelet_transform(chrom_depth, bins, chromosomes)
rpca <- robust_pca(wt$mat.wavelet.transform)

residuals <- inv_wavelet_transform(
  rpca$Sparse_Signal,
  wt$chrom.informed.wavelet
)
expected <- inv_wavelet_transform(
  rpca$Expected_Normal,
  wt$chrom.informed.wavelet
)

colnames(residuals) <- rownames(wt$mat.wavelet.transform)
colnames(expected) <- rownames(wt$mat.wavelet.transform)

segment_table <- segment_residuals(residuals, alpha = 0.005)
cn_wavelet <- assign_cn_state(chrom_depth, expected, segment_table)
```

### Chromosome-level negative-binomial workflow

```r
library(CHASM)

positions <- colnames(read_depth_chrom)[!(colnames(read_depth_chrom) %in% c("ID", "celltype"))]

cn_nb <- assign_cn_state.chrom(
  read_depth_chrom,
  positions,
  min_lib_size = 5000,
  max_lib_size = 60000,
  coverage_bin_size = 0.01,
  var_cut_off = 0.01
)
```

### Merge wavelet and negative-binomial calls

```r
library(CHASM)

cn_merged <- merge_calls(cn_wavelet, cn_nb)
```

The merged output includes `cn_state_final`, which combines segment-level and
chromosome-level evidence.

## Tutorial Scripts

Example scripts are available in [`tutorial/`](tutorial):

- `test_wavelet_subchrom.R`
- `test_negbinom_chrom.R`
- `run_chasm_full_pipeline.R`

These scripts are local usage examples and currently contain dataset-specific
paths that should be updated for a different environment.

An HTML tutorial version of the bundled full example pipeline is available at
[`tutorial/chasm-example-full-pipeline.html`](tutorial/chasm-example-full-pipeline.html),
with vignette source in
[`vignettes/chasm-example-full-pipeline.Rmd`](vignettes/chasm-example-full-pipeline.Rmd).
