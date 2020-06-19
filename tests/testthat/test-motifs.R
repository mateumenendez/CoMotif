context("Testing motifs")

test_that("motifs", {
  data(pos)
  data(neg)
  res <- motifs(pos = pos, neg = neg, network_name = "test_network", out = "count")
  groups <- res$groups
  expect_equal(sum(groups[,-1]), 42)
  count <- res$count
  expect_equal(sum(count[,2:8]), 14)
})
