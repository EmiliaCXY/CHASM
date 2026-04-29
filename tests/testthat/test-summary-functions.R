test_that("summarize_genome_altered computes altered percentage from chrom_bin calls", {
  cn_calls <- data.frame(
    ID = c("cell1", "cell1", "cell1", "cell2", "cell2", "cell2"),
    segment_id = c(
      "chr1_1_200", "chr1_1_200", "chr1_201_300",
      "chr1_1_200", "chr1_1_200", "chr1_201_300"
    ),
    chrom_bin = c(
      "chr1_1_100", "chr1_101_200", "chr1_201_300",
      "chr1_1_100", "chr1_101_200", "chr1_201_300"
    ),
    cn_state_final = c(3, 3, 2, 2, 2, 2),
    stringsAsFactors = FALSE
  )

  summary_df <- summarize_genome_altered(cn_calls)

  expect_equal(summary_df$genome_size, c(300, 300))
  expect_equal(summary_df$altered_length, c(200, 0))
  expect_equal(summary_df$pct_genome_altered, c(200 / 300 * 100, 0))
})

test_that("summarize_genome_altered accepts external bin_ids and cn_state_adj", {
  cn_calls <- data.frame(
    ID = c("cell1", "cell1", "cell2", "cell2"),
    segment_id = c("chr2_1_50", "chr2_51_100", "chr2_1_50", "chr2_51_100"),
    cn_state_adj = c(1, 2, 2, 2),
    stringsAsFactors = FALSE
  )

  summary_df <- summarize_genome_altered(
    cn_calls,
    bin_ids = c("chr2_1_50_p", "chr2_51_100_q")
  )

  expect_equal(summary_df$altered_length, c(50, 0))
  expect_equal(summary_df$genome_size, c(100, 100))
  expect_equal(summary_df$pct_genome_altered, c(50, 0))
})

test_that("prepare_cn_heatmap_input orders bins and inserts chromosome gaps", {
  cn_calls <- data.frame(
    ID = c("cell1", "cell1", "cell1", "cell2", "cell2", "cell2"),
    chrom_bin = c(
      "chr2_51_100", "chr1_1_50", "chr2_1_50",
      "chr2_51_100", "chr1_1_50", "chr2_1_50"
    ),
    cn_state_final = c(5, 2, 1, 4, 2, 0),
    stringsAsFactors = FALSE
  )

  heatmap_input <- .prepare_cn_heatmap_input(
    cn_calls,
    chromosomes = c("chr1", "chr2"),
    max_cn = 4,
    min_cn = 0
  )

  expect_equal(colnames(heatmap_input$matrix), c("chr1_1_50", "chr2_1_50", "chr2_51_100"))
  expect_equal(heatmap_input$gaps_col, 1)
  expect_equal(unname(heatmap_input$matrix["cell1", ]), c(2, 1, 4))
  expect_equal(unname(heatmap_input$matrix["cell2", ]), c(2, 0, 4))
})
