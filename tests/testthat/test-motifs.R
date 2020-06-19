context("Testing motifs")

test_that("motifs", {
  data(pos)
  data(neg)
  res <- motifs(pos = pos, neg = neg, network_name = "test_network", out = "count", square.motifs = T)
  groups <- res$groups
  expect_equal(sum(groups[,-1]), 102)
  count <- res$count
  expect_equal(sum(count[,2:8]), 14)
})
