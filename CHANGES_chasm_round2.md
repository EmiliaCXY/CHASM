# Changes on branch `chasm_round2`

Development log for the CHASM manuscript revision. Base: `main` @ `7af5239` ("update merge call criteria", 2026-05-06),
the version used for all revision results so far.

Goals of the branch:
1. Lock in the current default behaviour with a regression test.
2. Add the revision's methods as **opt-in** options (default behaviour must not change):
   - V2 back-transform for copy-number assignment;
   - second round: candidate regions, binomial seed support, masked SVD + regression on confident negatives
     (cells with no Strict CHASM call), hybrid background shift on supported candidates, strict re-segmentation.

Conventions:
- Every change to `R/` must keep `tests/testthat/test-regression-pipeline.R` passing (default pipeline = 7af5239 calls).
- Tests run with the `chasm_bioc` R (4.4.3) and testthat from a separate library `../Rlib`, with
  `NOT_CRAN=true` and BLAS threads set to 1, in a 1-CPU Slurm job (not on the login node).

---

## 2026-11-26

### Added
- `tests/testthat/helper-pipeline.R`
  - `run_reference_pipeline()`: the package's default pipeline, as in `tutorial/run_chasm_example_full_pipeline.R`
    (normalize → wavelet → robust PCA → inverse wavelet → `segment_residuals(alpha = 0.005)` → `assign_cn_state`
    → `assign_cn_state.chrom` with celltype "unknown" → `merge_calls`). It runs on a fixed subset of 300 cells
    (first 300 IDs, radix-sorted so the subset does not depend on locale) of the bundled `Read_depth_matrix_chr1_50mb_2_5pct_Rep1.rds`,
    with `set.seed(1234)` before robust PCA. It takes about 0.5 min on 1 CPU.
  - It runs the pipeline in a **fresh R process** (`callr::r`) with `LC_COLLATE=en_US.UTF-8`, sourcing `R/` from the package root,
    because CHASM's output depends on the collation set-up (see Found below).
  - NB-track cell IDs are normalised with `normalize_dnacopy_ids()`, as in the wavelet track. Without this, `merge_calls()` matches
    no cells and returns an empty table (found while building the fixture).
  - `compact_outputs()` keeps only the columns that define the calls:
    - `cn_bin`: ID, variable, segment_id, indicator, cn_state_raw, cn_state_adj;
    - `nb`: ID, chrom, p_value_adj, called_cna, cn_state_binom;
    - `final`: ID, variable, segment_id, cn_state_final.

    Observed and expected counts are dropped to keep the fixture small: 0.13 MB, against 3.5 MB with them.
  - `calls_key()`, `chasm_pkg_root()`, `regression_fixture_path()`: helpers.
- `tests/testthat/fixtures/make_regression_fixture.R`: one-time script that writes the reference output
  `tests/testthat/fixtures/regression_default_pipeline.rds` (xz-compressed). It records the git commit, whether `R/` was
  modified, R / DNAcopy / wavethresh versions, cell count, seed, alpha and runtime.
  The fixture was generated from the unmodified `R/` at 7af5239 with R 4.4.3, DNAcopy 1.80.0 and wavethresh 4.7.3.
  **Do not regenerate it after changing `R/`** unless a change to the default output is intended and reviewed.
- `tests/testthat/test-regression-pipeline.R`: reruns the reference pipeline and checks that:
  - columns and row counts match the reference;
  - merged calls are non-empty;
  - all kept columns are equal (tolerance 1e-8, order-independent).

  Skipped on CRAN.
- `CHANGES_chasm_round2.md` (this file).

### Test result
- With the fixture from unmodified 7af5239, the regression test passes (10 expectations; 1-CPU Slurm job, 2026-11-26).
  The only failures are the 5 pre-existing ones below.

### Found (not changed)
- **CHASM's output depends on string collation.**
  - **Cause:** `normalize_depth()` ends with `merge(..., by = "bin")`, whose row order follows the collation. That order sets the
    cell (column) order of the wavelet matrix, and the randomized robust PCA (`rsvd::rrpca`) gives a different background for a
    different cell order, even with the same seed.
  - **Size:** on the 300-cell example, the run under testthat's `LC_COLLATE=C` differs from the run under `en_US.UTF-8` in
    the robust-PCA background (max |difference| 0.027 in L), raw CN, and 3 cell × bin rows of the segmentation.
  - **Testing note:** resetting the locale with `Sys.setlocale()` inside a test did not restore the outside behaviour (ICU collation),
    hence the separate-process helper.
  - **Revision results are unaffected among themselves:** they all ran under `en_US.UTF-8`.
  - **Possible fix, not done:** make the order explicit, i.e. sort cells and bins by radix or by position before the wavelet
    transform. That would change default output under at least one locale, so it needs a reviewed fixture update.
- The existing test suite already fails on unmodified 7af5239 (5 problems):
  - `test-merge-calls.R:21`: `merge_calls()` errors ("replacement has 0 rows") on the test's toy input.
  - `test-negbinom.R:52` and `:69`: expected output columns of `cn_test_nb()` / `assign_cn_state.chrom()` not all present.
  - `test-summary-functions.R:71` and `:77`: expect the example file `1-27Nx_peakfrag_cov_5mb.rds`, which is not bundled
    in `inst/extdata`.

  Left as is for now; to be fixed separately.

### Not changed
- No changes to `R/`, `NAMESPACE`, `DESCRIPTION` or documentation yet.
