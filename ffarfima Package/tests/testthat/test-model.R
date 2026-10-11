sim <- sim_ffarfima(n = 160, d = 0.3, p = 1, q = 0, no_grid = 21, seed = 1)
X <- sim$X

test_that("every memory estimator returns a finite value", {
  for (m in c("lw", "dfa", "lpwn1", "lpwn2")) {
    expect_true(is.finite(estimate_memory(X, method = m)$d))
  }
  bc <- estimate_memory(X, method = "bc_lpwn2", B = 3, burn_in = 20, seed = 2)
  expect_equal(bc$d, bc$d_initial - bc$bias)
  expect_equal(bc$d, estimate_memory(X, method = "bc_lpwn2", B = 3, burn_in = 20, seed = 2)$d)
})

test_that("order selection stays inside the candidate set", {
  for (m in c("ic_full", "ic_auto")) {
    o <- select_order(X, d = 0.3, method = m, p = 1:2, q = 0:1, trim = 20)
    expect_true(o$p %in% 1:2 && o$q %in% 0:1)
    expect_true(is.integer(o$p) && is.integer(o$q))
  }
  cv <- select_order(X[, 1:100], d = 0.3, method = "cv", p = 0:1, q = 0:1, trim = 20,
                     valid = X[, 101:130], h = 3)
  expect_equal(nrow(cv$table), 3)
})

test_that("predict returns bands and an optional test evaluation", {
  fit <- ffarfima(X[, 1:100], X[, 101:130], d = 0.3, order = c(1, 0), h = 3, trim = 20, verbose = FALSE)
  fc <- predict(fit, verbose = FALSE)
  expect_equal(dim(fc$point_forecast), c(21, 3))
  expect_null(fc$evaluation)
  b <- fc$bands$SR[["80%"]]
  expect_true(all(b$lower <= fc$point_forecast & fc$point_forecast <= b$upper))
  fc2 <- predict(fit, test = X[, 131:160], verbose = FALSE)
  expect_equal(nrow(fc2$evaluation$point), 3)
  expect_setequal(unique(fc2$evaluation$interval$band), c("SR", "NS", "AC"))
  expect_error(predict(fit, h = 4, verbose = FALSE), "calibrated for h <= 3")
  expect_error(predict(fit, level = 0.9), "set in ffarfima")
  expect_identical(fit$order, c(p = 1L, q = 0L))
  expect_equal(dim(residuals(fit)), c(21, 80))
  expect_true(all(is.na(residuals(fit)[, 1])))
  fit_ac <- ffarfima(X[, 1:100], X[, 101:130], d = 0.3, order = c(1, 0), h = 3, band = "AC", trim = 20, verbose = FALSE)
  expect_equal(names(predict(fit_ac, verbose = FALSE)$bands), "AC")
  expect_equal((fitted(fit) + residuals(fit))[, -1], X[, 21:100][, -1], ignore_attr = TRUE)
  expect_equal(dim(fitted(fit, scale = "differenced")), c(21, 80))
  expect_error(predict(fit_ac, band = "SR", verbose = FALSE), "not chosen in ffarfima")
  expect_equal(nrow(as.data.frame(fc2, test = TRUE)), (28 * 3 + 2 + 1) * 21 * 3 * 2)
  grDevices::pdf(NULL)
  expect_silent(plot(fc2, h = 2, target_test_curve = 5))
  expect_silent(plot(fc2, h = 3))
  grDevices::dev.off()
  expect_error(plot(fc2, h = 3, target_test_curve = 2), "target test curves 3 to 30")
  r <- Filter(function(z) length(z$target_test_curve) >= 2 && z$target_test_curve[2] == 5, attr(fc2, "test_forecasts"))[[1]]
  expect_named(fc2$evaluation, c("point", "interval", "wis"))
  expect_named(fc2, c("point_forecast", "bands", "evaluation", "origin", "settings", "model"))
  expect_equal(r$observed[, 2], X[, 130 + 5])
  bench <- ffarfima(X[, 1:100], X[, 101:130], d = 1, method = "far", h = 3, trim = 20, verbose = FALSE)
  fb <- predict(bench, test = X[, 131:160], verbose = FALSE)
  expect_true(all(is.finite(fb$evaluation$point$RMSPE)))
})

test_that("blocks out of time order give a warning", {
  Y <- X
  colnames(Y) <- 1850 + seq_len(ncol(Y))
  fit <- ffarfima(Y[, 1:100], d = 0.3, order = c(1, 0), trim = 20, verbose = FALSE)
  expect_warning(predict(fit, test = Y[, 131:160], h = 3, verbose = FALSE), "must follow each other")
  expect_silent(predict(fit, test = Y[, 101:130], h = 3, verbose = FALSE))
})

test_that("a fit without validation gives point forecasts only", {
  fit <- ffarfima(X[, 1:100], d = 0.3, order = c(1, 0), trim = 20, verbose = FALSE)
  fc <- predict(fit, test = X[, 101:130], h = 3, verbose = FALSE)
  expect_length(fc$bands, 0)
  expect_null(fc$evaluation$interval)
  expect_equal(nrow(fc$evaluation$point), 3)
})
