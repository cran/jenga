test_that("multiscale settings are normalized", {
  expect_equal(jenga:::normalize_multiscale(c(4, 2, 2, 1)), c(1L, 2L, 4L))
  expect_equal(jenga:::normalize_multiscale(c(NA, 0, -1)), 1L)
})

test_that("multiscale embeddings preserve rows and add aggregate blocks", {
  x <- data.frame(a = 1:8, b = 9:16, c = 17:24, d = 25:32)
  embedded <- jenga:::multiscale_embed(x, c(1, 2, 4))

  expect_equal(nrow(embedded), nrow(x))
  expect_gt(ncol(embedded), ncol(x))
  expect_false(any(!is.finite(embedded)))
  expect_true(any(grepl("^s2_", colnames(embedded))))
  expect_true(any(grepl("^s4_", colnames(embedded))))
})

test_that("projection ANN returns stable candidates without changing RNG state", {
  embedding <- matrix(seq_len(60), nrow = 12)

  set.seed(123)
  old_seed <- .Random.seed
  candidates <- jenga:::projection_candidates(
    jenga:::normalize_embedding(embedding),
    candidate_count = 5,
    reference_id = 12,
    seed = 42
  )

  expect_identical(.Random.seed, old_seed)
  expect_length(candidates, 5)
  expect_true(12 %in% candidates)
  expect_true(all(candidates >= 1 & candidates <= nrow(embedding)))
})

test_that("nearest-neighbor candidate backends honor exact and projection modes", {
  sequences <- list(
    data.frame(t1 = 1:40, t2 = 2:41, t3 = 3:42),
    data.frame(t1 = 40:1, t2 = 41:2, t3 = 42:3)
  )

  exact <- jenga:::neighbour_candidates(
    sequences,
    method = c("euclidean", "euclidean"),
    k = 3,
    ann = "exact",
    multiscale = c(1, 2),
    seed = 1,
    reference_id = 40
  )
  projection <- jenga:::neighbour_candidates(
    sequences,
    method = c("euclidean", "euclidean"),
    k = 3,
    ann = "projection",
    multiscale = c(1, 2),
    seed = 1,
    reference_id = 40
  )

  expect_equal(exact, seq_len(40))
  expect_true(40 %in% projection)
  expect_true(length(projection) < length(exact))
  expect_true(all(projection %in% exact))
})

test_that("conformal prediction adds calibrated interval columns", {
  quant_pred <- data.frame("50%" = c(10, 20, 30), check.names = FALSE)
  raw_errors <- matrix(c(1, 2, 3, 2, 3, 4, 1, 4, 5), nrow = 3)

  calibrated <- jenga:::add_conformal_prediction(
    quant_pred,
    raw_errors = raw_errors,
    ci = 0.8,
    binary_class = FALSE
  )

  expect_named(
    calibrated,
    c("50%", "conformal_lower", "conformal_upper", "conformal_radius", "conformal_coverage")
  )
  expect_equal(calibrated$conformal_coverage, rep(0.8, 3))
  expect_true(all(calibrated$conformal_lower <= quant_pred[["50%"]]))
  expect_true(all(calibrated$conformal_upper >= quant_pred[["50%"]]))
})

test_that("binary conformal intervals are clipped to class bounds", {
  quant_pred <- data.frame("50%" = c(0.2, 0.8), check.names = FALSE)
  raw_errors <- matrix(c(0.5, 0.7, 0.5, 0.7), nrow = 2)

  calibrated <- jenga:::add_conformal_prediction(
    quant_pred,
    raw_errors = raw_errors,
    ci = 0.9,
    binary_class = TRUE
  )

  expect_true(all(calibrated$conformal_lower >= 0))
  expect_true(all(calibrated$conformal_upper <= 1))
})

test_that("positive-only histories keep conformal and plot bands nonnegative", {
  positive_history <- c(4, 5, 6, 7)
  quant_pred <- data.frame(
    min = c(-3, 1),
    "10%" = c(-1, 1.5),
    "50%" = c(2, 3),
    "90%" = c(4, 5),
    max = c(5, 6),
    mean = c(2.5, 3.5),
    mode = c(-0.5, 3),
    check.names = FALSE
  )
  raw_errors <- matrix(c(10, 10, 10, 10), nrow = 2)

  calibrated <- jenga:::add_conformal_prediction(
    quant_pred,
    raw_errors = raw_errors,
    ci = 0.8,
    binary_class = FALSE,
    ts = positive_history
  )
  plot_data <- jenga:::plotter(
    calibrated,
    ci = 0.8,
    ts = positive_history,
    dates = NULL,
    feat_name = "positive_series"
  )

  expect_true(all(calibrated$conformal_lower >= 0))
  expect_true(all(calibrated[["min"]] >= 0))
  expect_true(all(calibrated[["10%"]] >= 0))
  expect_true(all(calibrated[["mode"]] >= 0))
  expect_true(all(plot_data$lower >= 0))
  expect_true(all(plot_data$y_forcat >= 0))
})

test_that("binary distance metrics are finite and named", {
  metrics <- jenga:::binary_distance_metrics(
    holdout = c(1, 0, 1, 1, NA),
    forecast = c(1, 1, 0, 1, 1)
  )

  expect_named(
    metrics,
    c("dice", "jaccard", "cosine", "canberra", "gower", "tanimoto", "hassebrook", "taneja", "lorentzian", "clark", "sorensen", "harmonic_mean", "avg")
  )
  expect_true(all(is.finite(metrics)))
  expect_equal(jenga:::binary_distance_metrics(c(1, 0), c(1, 0))[["dice"]], 0)
  expect_equal(jenga:::binary_distance_metrics(c(0, 0), c(0, 0))[["jaccard"]], 0)
})

test_that("numeric error metrics preserve the expected metric contract", {
  metrics <- jenga:::numeric_error_metrics(
    holdout = c(1, 2, 3),
    forecast = c(1, 2, 4),
    benchmark = c(1, 1, 1),
    scale = 1
  )

  expect_named(
    metrics,
    c("me", "mae", "mse", "rmsse", "mpe", "mape", "rmae", "rrmse", "rame", "mase", "smse", "sce", "gmrae")
  )
  expect_true(all(is.finite(metrics[names(metrics) != "gmrae"])))
})

test_that("plotter returns dependency-free jenga plot objects", {
  plot_data <- jenga:::ts_graph(
    x_hist = 1:5,
    y_hist = 1:5,
    x_forcat = 6:8,
    y_forcat = c(6, 7, 8),
    lower = c(5, 6, 7),
    upper = c(7, 8, 9)
  )

  expect_s3_class(plot_data, "jenga_plot")
})

test_that("expanded distance methods are available and finite", {
  methods <- c("euclidean", "manhattan", "minkowski", "chebyshev", "canberra", "clark", "sorensen", "lorentzian", "cosine", "correlation", "hamming")
  sequence <- data.frame(t1 = c(1, 2, 3, 4), t2 = c(2, 2, 5, 4), t3 = c(3, 1, 6, 4))

  for(method in methods)
  {
    distances <- jenga:::sequence_distances(sequence, method, candidate_id = 1:4, reference_id = 4, gpu = "off")
    expect_length(distances, 4)
    expect_true(all(is.finite(distances)))
    expect_equal(distances[4], 0)
  }
})

test_that("expanded kernels are validated and produce finite forecasts", {
  expect_equal(jenga:::supported_kernels(), c("norm", "cauchy", "logis", "unif", "t", "laplace", "epanechnikov", "triangular", "biweight", "triweight", "cosine"))
  expect_equal(jenga:::supported_distance_methods(), c("euclidean", "manhattan", "minkowski", "chebyshev", "canberra", "clark", "sorensen", "lorentzian", "cosine", "correlation", "hamming"))

  expect_error(jenga:::validate_choices("bad_kernel", jenga:::supported_kernels(), "kernel"), "kernel must be one of")
  expect_equal(jenga:::kernel_weighted_mean(c(1, 2, 3), c(0, 0, 0)), 2)
})

test_that("reconstructed sign filter removes invalid signs before quantiles", {
  positive_history <- c(1.2, 3.1, 2.4, 4.5)
  reconstructed <- matrix(
    c(-10, -5, 1, 2, -3, -2, -1, -4),
    nrow = 4,
    ncol = 2
  )

  filtered <- jenga:::doxa_filter_plus(
    positive_history,
    reconstructed,
    binary_class = FALSE,
    seed = 42,
    sign_filter = "remove"
  )

  expect_true(any(is.na(filtered[, 1])))
  expect_true(all(filtered[, 1][is.finite(filtered[, 1])] >= 0))
  expect_true(all(filtered[, 2] == 0))

  quantiles <- jenga:::fast_qpred(
    raw_pred = filtered,
    ts = positive_history,
    ci = 0.8,
    binary_class = FALSE,
    dates = NULL,
    seed = 42
  )

  expect_true(all(quantiles[["min"]] >= 0))
})

test_that("reconstructed sign filter respects negative-only histories", {
  negative_history <- c(-5, -3, -4)
  reconstructed <- matrix(c(-6, 2, -1, 3), nrow = 2)

  filtered <- jenga:::doxa_filter_plus(
    negative_history,
    reconstructed,
    binary_class = FALSE,
    seed = 42,
    sign_filter = "remove"
  )

  expect_true(any(is.na(filtered)))
  expect_true(all(filtered[is.finite(filtered)] <= 0))
})
