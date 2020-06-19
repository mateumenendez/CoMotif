context("Testing conf.motifs")

test_that("conf.motifs", {
  data(pos)
  data(neg)
  RNGversion("3.6")
  set.seed(999)
  res <- conf.motifs(pos = pos, neg = neg, network_name = "test_network",
                     num_random_networks = 10, out = "normalized")
  expect_equal(round(sum(res[,2:8]), digits = 2), 1.85)
})
