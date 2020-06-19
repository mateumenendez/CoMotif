context("Testing conf.motifs.p")

test_that("conf.motifs.p", {
  data(pos)
  data(neg)
  RNGversion("3.6")
  set.seed(999)
  res <- conf.motifs.p(pos = pos, neg = neg, network_name = "test_network",
                     num_random_networks = 12, out = "count", cores = 2, square.motifs = T)
  expect_equal(round(sum(res[,14:15]), digits = 2), 408)
})
