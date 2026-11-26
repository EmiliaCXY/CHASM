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

## 2026-11-26 (3): V2 back-transform option in `assign_cn_state()`

### Added
- `R/cn_assign.R`
  - `assign_cn_state(..., backtransform = c("V0", "V2"))`, default `"V0"` (unchanged output).
  - New internal `backtransform_expected_v2()`:
    `e_V2 = lib × ((√(e_V0/lib) − β̂_med + β̂_mean)² + Var_b)`.
    - β̂_med, β̂_mean and Var_b are the per-bin median, mean and population variance of `beta_i_c` over the reference cells
      (library ≥ 75th percentile), the same cells that define `beta_i_hat`.
    - **Why:** V0 squares a median of square-root shares, which underestimates sparse bins and shallow cells, since E[x] = (E√x)² + Var(√x).
      That inflates observed/expected and gives false gains, e.g. on chrX and small autosomes.
- `tests/testthat/test-backtransform-v2.R`:
  1. With identical segmentation, V2 expected counts, raw CN and adjusted CN equal an independent re-implementation of the
     revision formula (`call_round()`; tolerance 1e-10).
  2. On chromosomes the example spike-in does not touch, the median raw CN of V2 is closer to 2 than V0's.
- `tests/testthat/helper-pipeline.R`: `run_reference_pipeline()` gains `backtransform` and `compact` arguments.
  The default path is unchanged.

### Test result (1-CPU Slurm job)
- The V2 tests pass (7 expectations). The regression test still passes, so the default calls are identical to 7af5239.
  The same 5 pre-existing failures remain.

### Evaluation: does V2 behave better than V0? (`revision_claude/Scripts/task8_chasm_package/01_eval_v2_package.R`)
- **Set-up:** this branch's package code was applied to the saved CHASM runs of 160 spike-in settings, Rep1:
  - breast 5/10/20/50 Mb and chr1p, kidney 5/10/20 Mb;
  - CN ±1/±2;
  - 1/5/10/15/20% prevalence.

  Robust PCA and segmentation are reused, so only CN assignment differs. Package V0 reproduces the saved calls in all 160 settings (100%).
- Medians over settings, V0 → V2:

  | | Breast, before NB merge | Breast, final calls | Kidney, before NB merge | Kidney, final calls |
  |---|---|---|---|---|
  | Bin FPR | 9.9% → 1.4% | 0.45% → 0.34% | 1.26% → 0.65% | 0.50% → 0.48% |
  | Normal cells with any call | 84% → 42% | 25.6% → 25.3% | 39% → 31% | 29.7% → 28.9% |
  | Normal cells with > 5% of genome called | 74% → 13% | 2.4% → 1.2% | 12.5% → 4.5% | 2.0% → 1.9% |
  | Gains on chrX/13/18/21 per normal cell | 1.26 → 0.22 | 0.058 → 0.045 | 0.21 → 0.08 | 0.057 → 0.041 |
  | Carriers detected | 0.425 → 0.417 | 0.407 → 0.412 | 0.204 → 0.198 | 0.204 → 0.198 |

- **Bias** (median raw CN of normal cells' single-segment chromosomes, V0 → V2):
  - breast: 2.18 → 2.07 overall, chrX 2.66 → 2.07;
  - kidney: 2.08 → 2.02 overall, chrX 2.25 → 2.07.
- **Final calls, all 160 settings:**
  - V2 lowers the bin FPR in 160 of 160 settings, and lowers the share of normal cells with > 5% of the genome called in 144.
  - Detection is within ±0.02 in 143. It is lower in 13, all of them gains: worst kidney 20 Mb +2 20% 0.57 → 0.47 and breast chr1p +1 20% 0.28 → 0.20.
    It is higher in 4, all of them one-copy losses: e.g. breast chr1p −1 20% 0.51 → 0.59.
- **Trade-off:** false losses on sparse chromosomes rise slightly (final, median +0.011 per normal cell), because V2 over-corrects shallow cells
  (the depth dependence noted in the revision).
- **Verdict:** V2 removes most of V0's sparse-bin gain bias and cuts genome-wide false calls, at a small cost in detection for
  high-prevalence gains.

## 2026-11-26 (2): faster long-to-wide pivot in `wavelet_transform()`

### Changed
- `R/wavelet.R`, `wavelet_transform()`: the long → wide pivot (cells × bins) now uses `data.table::dcast()` instead of
  `stats::reshape(direction = "wide")`. It comes from `analysis_claude/Tools/chasm_patches.R` (`wavelet_transform_fast`),
  with one correction:
  - `dcast` sorts cells by name, while `reshape` keeps them in order of first appearance. The cell order is the column order of
    the wavelet matrix and changes the randomized robust PCA, so the patch as written would **not** have reproduced the calls.
  - Rows are therefore put back in first-appearance order (`match(unique(barcode), ...)`). Column names stay
    `read_depth_sqrt_centered.<bin>`, and the bin reordering by chromosome / arm / position is unchanged.
- `DESCRIPTION`: `data.table` added to Imports.
- The patch's other part, `robust_pca_fast()` (rsvd `maxiter` 20 instead of 50), is **not** included: it changes the output.

### Added
- `tests/testthat/test-wavelet-dcast.R`: on 120 cells with the long table shuffled, the new pivot gives the same cell order as
  the former reshape code, and the inverse wavelet transform reproduces the reshape matrix (tolerance 1e-12).

### Test result (1-CPU Slurm job)
- Regression test passes: the default-pipeline calls are identical to the 7af5239 reference. The new dcast test passes.
  The same 5 pre-existing failures remain.
- Speed (pivot step; old `reshape` alone vs the whole new `wavelet_transform()`), 1 CPU:

  | Cells | Old | New |
  |---|---|---|
  | 2,000 | 44.1 s | 1.2 s |
  | 10,000 | 234.7 s | 5.2 s |

  About 40× faster, and the cell order is identical at both sizes.

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
