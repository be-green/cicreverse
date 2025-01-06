
testthat::test_that("Counterfactuals are correct", {
  test_data = data.table::data.table(
    y = c(1, 0, 0, 0),
    t = c(1, 1, 0, 0),
    g = c(1, 1, 1, 1),
    i = c(1, 2, 1, 2)
  )
  cfs = get_all_cfs(i = test_data$i,t = test_data$t,
                    y = test_data$y, g = test_data$g)

  testthat::expect_equal(unique(cfs$i), 1)
  testthat::expect_equal(unique(setdiff(cfs$j, cfs$i)), c(2))
  testthat::expect_equal(unique(sort(cfs$s)), c(0, 1))
})

testthat::test_that("Output contains correct counterfactuals", {

  set.seed(42)

  n = 100
  y = c(rep(1, n), rep(0, 5 * n))
  t = c(rep(1, n * 3), rep(0, n * 3))
  g = rep(1, n * 6)
  i = c(rep(1, n), rep(0, n), rep(2, n), rep(1, n), rep(0, n), rep(2, n))
  x = rnorm(n * 6)
  tmp = cic(i, t, y, x, g, t_window = NA, n_boots = 20)

  testthat::expect_equal(unique(tmp$cic_data$j), c(0, 2))
  testthat::expect_equal(unique(tmp$cic_data$t), 1)
  testthat::expect_equal(unique(tmp$cic_data$s), c(0))
  testthat::expect_equal(unique(tmp$training_data$i), 1)
  testthat::expect_equal(unique(tmp$training_data$t), 1)
  testthat::expect_equal(max(tmp$bootstraps$b), 20)
})

# plot_data = combined_cf(tmp)
