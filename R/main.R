###
as_mapper <- function(.f)
{
  if(inherits(.f, "formula"))
  {
    f_env <- environment(.f)
    f_expr <- .f[[length(.f)]]
    return(function(...){
      args <- list(...)
      eval_env <- new.env(parent = f_env)
      if(length(args) >= 1)
      {
        assign(".x", args[[1]], envir = eval_env)
        assign(".", args[[1]], envir = eval_env)
      }
      if(length(args) >= 2){assign(".y", args[[2]], envir = eval_env)}
      for(i in seq_along(args)){assign(paste0("..", i), args[[i]], envir = eval_env)}
      eval(f_expr, envir = eval_env)
    })
  }
  return(.f)
}

###
map <- function(.x, .f, ...)
{
  .f <- as_mapper(.f)
  lapply(.x, .f, ...)
}

###
map_lgl <- function(.x, .f, ...)
{
  vapply(map(.x, .f, ...), isTRUE, logical(1))
}

###
map_dbl <- function(.x, .f, ...)
{
  vapply(map(.x, .f, ...), function(x) as.numeric(x)[1], numeric(1))
}

###
map_dfr <- function(.x, .f, ...)
{
  rows <- map(.x, .f, ...)
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if(length(rows) == 0){return(data.frame())}
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  as.data.frame(out)
}

###
map_dfc <- function(.x, .f, ...)
{
  cols <- map(.x, .f, ...)
  cols <- lapply(cols, as.data.frame)
  if(length(cols) == 0){return(data.frame())}
  as.data.frame(do.call(cbind, cols))
}

###
map2 <- function(.x, .y, .f, ...)
{
  .f <- as_mapper(.f)
  n <- min(length(.x), length(.y))
  lapply(seq_len(n), function(i) .f(.x[[i]], .y[[i]], ...))
}

###
map2_dbl <- function(.x, .y, .f, ...)
{
  vapply(map2(.x, .y, .f, ...), function(x) as.numeric(x)[1], numeric(1))
}

###
map2_dfc <- function(.x, .y, .f, ...)
{
  cols <- map2(.x, .y, .f, ...)
  cols <- lapply(cols, as.data.frame)
  if(length(cols) == 0){return(data.frame())}
  as.data.frame(do.call(cbind, cols))
}

###
pmap <- function(.l, .f, ...)
{
  .f <- as_mapper(.f)
  n <- min(vapply(.l, length, integer(1)))
  lapply(seq_len(n), function(i){
    args <- lapply(.l, function(x) x[[i]])
    names(args) <- NULL
    do.call(.f, c(args, list(...)))
  })
}

###
pmap_dbl <- function(.l, .f, ...)
{
  vapply(pmap(.l, .f, ...), function(x) as.numeric(x)[1], numeric(1))
}

###
transpose <- function(.x)
{
  if(length(.x) == 0){return(list())}
  n <- max(vapply(.x, length, integer(1)))
  lapply(seq_len(n), function(i) lapply(.x, function(x) x[[i]]))
}

###
split <- function(x, ..., along = NULL)
{
  if(is.null(along)){return(base::split(x, ...))}
  if(along == 1){return(lapply(seq_len(nrow(x)), function(i) x[i,, drop = FALSE]))}
  if(along == 2){return(lapply(seq_len(ncol(x)), function(i) x[, i, drop = FALSE]))}
  stop("only along = 1 or along = 2 is supported")
}

###
format_elapsed_time <- function(seconds)
{
  seconds <- max(0, round(seconds))
  hours <- seconds %/% 3600
  minutes <- (seconds %% 3600) %/% 60
  secs <- seconds %% 60
  sprintf("%02d:%02d:%02d", hours, minutes, secs)
}

###
new_jenga_plot <- function(x_hist, y_hist, x_forcat, y_forcat, lower = NULL, upper = NULL, label_x = "Horizon", label_y = "Forecasted Var", date_format = "%b-%d-%Y")
{
  structure(list(x_hist = x_hist, y_hist = y_hist, x_forcat = x_forcat, y_forcat = y_forcat, lower = lower, upper = upper, label_x = label_x, label_y = label_y, date_format = date_format), class = "jenga_plot")
}

#' Plot jenga forecast objects
#'
#' @param x A jenga forecast plot object.
#' @param ... Additional arguments passed to the base plotting device.
#'
#' @export
plot.jenga_plot <- function(x, ...)
{
  bounds <- history_level_bounds(x$y_hist)
  y_forcat <- apply_level_bounds(x$y_forcat, bounds)
  lower <- apply_level_bounds(x$lower, bounds)
  upper <- apply_level_bounds(x$upper, bounds)

  y_all <- c(as.numeric(x$y_hist), y_forcat, lower, upper)
  y_range <- range(y_all[is.finite(y_all)], na.rm = TRUE)
  if(!all(is.finite(y_range))){y_range <- c(0, 1)}
  graphics::plot(c(x$x_hist, x$x_forcat), c(as.numeric(x$y_hist), rep(NA_real_, length(x$x_forcat))), type = "l", col = "gray43", lwd = 1.3, xlab = x$label_x, ylab = x$label_y, ylim = y_range, ...)
  if(!is.null(lower) && !is.null(upper))
  {
    graphics::polygon(c(x$x_forcat, rev(x$x_forcat)), c(lower, rev(upper)), col = grDevices::adjustcolor("seagreen2", alpha.f = 0.3), border = NA)
  }
  graphics::lines(x$x_forcat, y_forcat, col = "seagreen4", lwd = 1.3)
  invisible(x)
}

#' Print jenga forecast plot objects
#'
#' @param x A jenga forecast plot object.
#' @param ... Additional arguments passed to the base plotting device.
#'
#' @export
print.jenga_plot <- function(x, ...)
{
  plot(x, ...)
  invisible(x)
}

###
is_date <- function(x)
{
  inherits(x, "Date")
}

###
dummy_cols_fast <- function(df)
{
  out <- list()
  for(col_name in names(df))
  {
    x <- df[[col_name]]
    x <- as.character(x)
    levels <- sort(unique(x[!is.na(x)]))
    if(length(levels) == 0){next}
    most_frequent <- names(which.max(table(x)))
    levels <- setdiff(levels, most_frequent)
    for(level in levels)
    {
      safe_name <- make.names(paste(col_name, level, sep = "_"), unique = TRUE)
      out[[safe_name]] <- as.integer(!is.na(x) & x == level)
    }
  }
  as.data.frame(out, check.names = FALSE)
}

###
fast_time_impute <- function(df)
{
  impute_one <- function(x)
  {
    x <- as.numeric(x)
    ok <- is.finite(x)
    if(all(ok)){return(x)}
    if(!any(ok)){return(rep(0, length(x)))}
    stats::approx(which(ok), x[ok], xout = seq_along(x), rule = 2, ties = "ordered")$y
  }
  as.data.frame(lapply(df, impute_one))
}

###
smooth_series <- function(x)
{
  x <- as.numeric(x)
  if(length(x) < 5 || length(unique(x[is.finite(x)])) < 3){return(x)}
  idx <- seq_along(x)
  span <- min(1, max(0.25, 10/length(x)))
  fitted <- tryCatch(stats::loess(x ~ idx, span = span, degree = 1, na.action = stats::na.exclude)$fitted, error = function(e) x)
  as.numeric(fitted)
}

###
sample_mode <- function(x)
{
  x <- x[is.finite(x)]
  if(length(x) == 0){return(NA_real_)}
  if(length(unique(x)) == 1){return(x[1])}
  bins <- max(10, ceiling(sqrt(length(x))))
  hist_data <- graphics::hist(x, breaks = bins, plot = FALSE)
  hist_data$mids[which.max(hist_data$counts)]
}

###
sample_skewness <- function(x)
{
  x <- x[is.finite(x)]
  if(length(x) < 3){return(NA_real_)}
  centered <- x - mean(x)
  s <- sd(x)
  if(!is.finite(s) || s == 0){return(0)}
  mean((centered/s)^3)
}

###
sample_kurtosis <- function(x)
{
  x <- x[is.finite(x)]
  if(length(x) < 4){return(NA_real_)}
  centered <- x - mean(x)
  s <- sd(x)
  if(!is.finite(s) || s == 0){return(0)}
  mean((centered/s)^4)
}

###
sample_entropy <- function(x)
{
  x <- x[is.finite(x)]
  if(length(x) == 0){return(NA_real_)}
  p <- as.numeric(table(x))/length(x)
  -sum(p * log(p))
}

###
supported_distance_methods <- function()
{
  c("euclidean", "manhattan", "minkowski", "chebyshev", "canberra", "clark", "sorensen", "lorentzian", "cosine", "correlation", "hamming")
}

###
supported_kernels <- function()
{
  c("norm", "cauchy", "logis", "unif", "t", "laplace", "epanechnikov", "triangular", "biweight", "triweight", "cosine")
}

###
history_level_bounds <- function(ts, binary_class = FALSE)
{
  if(binary_class){return(c(lower = 0, upper = 1))}

  ts <- as.numeric(ts)
  finite_history <- ts[is.finite(ts)]
  bounds <- c(lower = -Inf, upper = Inf)
  if(length(finite_history) == 0){return(bounds)}

  if(all(finite_history >= 0)){bounds["lower"] <- 0}
  if(all(finite_history <= 0)){bounds["upper"] <- 0}
  return(bounds)
}

###
apply_level_bounds <- function(x, bounds)
{
  if(is.null(x)){return(NULL)}
  x <- as.numeric(x)
  if(is.finite(bounds[["lower"]])){x <- pmax(bounds[["lower"]], x)}
  if(is.finite(bounds[["upper"]])){x <- pmin(bounds[["upper"]], x)}
  return(x)
}

###
prediction_level_columns <- function(quant_pred)
{
  percentile_columns <- grep("^[0-9.]+%$", colnames(quant_pred), value = TRUE)
  intersect(c("min", percentile_columns, "max", "mean", "mode", "prop", "conformal_lower", "conformal_upper"), colnames(quant_pred))
}

###
enforce_prediction_sign_bounds <- function(quant_pred, ts, binary_class = FALSE)
{
  bounds <- history_level_bounds(ts, binary_class)
  if(!any(is.finite(bounds))){return(quant_pred)}

  level_columns <- prediction_level_columns(quant_pred)
  if(length(level_columns) == 0){return(quant_pred)}

  quant_pred[level_columns] <- lapply(quant_pred[level_columns], apply_level_bounds, bounds = bounds)
  return(quant_pred)
}

###
validate_choices <- function(x, allowed, label)
{
  if(is.null(x)){return(invisible(TRUE))}
  bad <- setdiff(x, allowed)
  if(length(bad) > 0){stop(label, " must be one of: ", paste(allowed, collapse = ", "), call. = FALSE)}
  invisible(TRUE)
}

#' jenga
#'
#' @param df A data frame with time features on columns (numerical or categorical features, but not both).
#' @param seq_len Positive integer. Time-step number of the projected sequence
#' @param smoother Logical. Perform optimal smoothing using standard loess (only for numerical features). Default: FALSE
#' @param k Positive integer. Number of neighbors to consider when applying kernel average. Min number is 3. Default: NULL (automatic selection).
#' @param method String. Distance method for calculating neighbors. Possible options are: "euclidean", "manhattan", "minkowski", "chebyshev", "canberra", "clark", "sorensen", "lorentzian", "cosine", "correlation", and "hamming". Default: NULL (automatic selection).
#' @param kernel String. Kernel used to calculate neighbor weights. Possible options are: "norm", "cauchy", "logis", "unif", "t", "laplace", "epanechnikov", "triangular", "biweight", "triweight", and "cosine". Default: NULL (automatic selection).
#' @param ci Confidence interval. Default: 0.8
#' @param n_windows Positive integer. Number of validation tests to measure/sample error. Default: 10.
#' @param mode String. Sequencing method: deterministic ("segmented"), or non-deterministic ("sampled"). Default: NULL (automatic selection).
#' @param n_sample Positive integer. Number of samples for grid or random search. Default: 30.
#' @param search String. Two option available: "grid", "random". Default: "random".
#' @param dates Date. Vector with dates for time features.
#' @param error_scale String. Scale for the scaled error metrics. Two options: "naive" (average of naive one-step absolute error for the historical series) or "deviation" (standard error of the historical series). Default: "naive".
#' @param error_benchmark String. Benchmark for the relative error metrics. Two options: "naive" (sequential extension of last value) or "average" (mean value of true sequence). Default: "naive".
#' @param seed Positive integer. Random seed. Default: 42.
#' @param ann String. Approximate nearest-neighbor backend. Options are: "auto", "exact", and "projection". Default: "auto".
#' @param gpu String. Acceleration policy for distance calculations. Options are: "auto", "off", and "torch". Default: "auto".
#' @param probabilistic Logical. Return simulated forecast distributions and probabilistic summaries. Default: TRUE.
#' @param conformal Logical. Add split-conformal forecast intervals calibrated on validation windows. Default: TRUE.
#' @param multiscale Positive integer vector. Sequence aggregation scales used in the neighbor representation. Default: c(1, 2, 4).
#' @param n_draws Positive integer. Number of simulated draws for probabilistic forecasts. Default: 1000.
#'
#' @author Giancarlo Vercellino \email{giancarlo.vercellino@gmail.com}
#'
#' @return This function returns a list including:
#' \itemize{
#' \item exploration: list of all models, complete with predictions, test metrics, prediction stats and plot
#' \item history: a table with the sampled models, hyper-parameters, validation errors
#' \item best_model: results for the best model, including:
#' \itemize{
#' \item predictions: min, max, q25, q50, q75, quantiles at selected ci, and different statics for numerical and categorical variables
#' \item testing_errors: training and testing errors for one-step and sequence for each ts feature (different measures for numerical and categorical variables)
#' }
#' \item distribution: simulated forecast draws for each feature when probabilistic = TRUE
#' \item settings: v2 execution settings for ANN, GPU, conformal, probabilistic and multiscale options
#' \item time_log
#' }
#'
#' @export
#'
#' @importFrom stats approx ecdf lm loess na.exclude pf quantile rnorm sd
#' @importFrom utils head tail
#' @importFrom stats dcauchy dlogis dnorm dt dunif
#' @examples
#' jenga(covid_in_europe[, c(2, 3)], n_sample = 1)
#' jenga(covid_in_europe[, c(4, 5)], n_sample = 1)
jenga <- function(df, seq_len = NULL, smoother = FALSE, k = NULL, method = NULL, kernel = NULL, ci = 0.8, n_windows = 10, mode = NULL, n_sample = 30, search = "random", dates = NULL, error_scale = "naive", error_benchmark = "naive", seed = 42, ann = "auto", gpu = "auto", probabilistic = TRUE, conformal = TRUE, multiscale = c(1, 2, 4), n_draws = 1000)
{
  set.seed(seed)

  start_time <- proc.time()[["elapsed"]]

  ann <- match.arg(ann, c("auto", "exact", "projection"))
  gpu <- match.arg(gpu, c("auto", "off", "torch"))
  if(!is.logical(probabilistic) || length(probabilistic) != 1){stop("probabilistic must be TRUE or FALSE")}
  if(!is.logical(conformal) || length(conformal) != 1){stop("conformal must be TRUE or FALSE")}
  n_draws <- as.integer(n_draws[1])
  if(is.na(n_draws) || n_draws < 1){stop("n_draws must be a positive integer")}
  multiscale <- normalize_multiscale(multiscale)
  if(gpu == "torch" && !torch_gpu_ready()){message("torch GPU is not available; falling back to CPU distances\n")}

  if(!is.data.frame(df)){stop("time features must be in dataframe format")}

  n_ts <- nrow(df)

  class_index <- any(map_lgl(df, ~ is.factor(.x) | is.character(.x)))
  all_classes <- all(class_index)
  numeric_index <- map_lgl(df, ~ is.integer(.x) | is.numeric(.x))
  all_numerics <- all(numeric_index)
  if(!(all_classes | all_numerics)){stop("only all numerics or all classes, not both")}

  if(all_classes){df <- dummy_cols_fast(df); binary_class <- rep(T, ncol(df))}
  if(all_numerics){binary_class <- rep(FALSE, ncol(df))}

  if(anyNA(df) & all_numerics){df <- fast_time_impute(df); message("linear time imputation on time features\n")}
  if(anyNA(df) & all_classes){df <- floor(fast_time_impute(df)); message("linear time imputation on time features\n")}
  if(smoother == TRUE & all_numerics){df <- as.data.frame(map(df, smooth_series)); message("performing loess smoothing\n")}

  n_feat <- ncol(df)
  feat_names <- colnames(df)
  distance_methods <- supported_distance_methods()
  kernels <- supported_kernels()
  validate_choices(method, distance_methods, "method")
  validate_choices(kernel, kernels, "kernel")

  deriv <- map_dbl(df, ~ best_deriv(.x))
  fixed <- FALSE
  if(all(c(length(seq_len)==1, length(k)==1, length(method)==n_feat, length(kernel)==n_feat, length(mode)==1))){fixed <- TRUE}

  if(!is.null(k) && any(k < 3)){k[k < 3] <- 3; message("setting minimum k value to 3")}
  if(!is.null(seq_len) && any(seq_len < 3)){seq_len[seq_len < 3] <- 3; message("setting minimum seq_len value to 3")}

  ###SAMPLERS
  sl_range <- c(3:floor(sqrt(n_ts)))
  k_range <- c(3:floor(sqrt(n_ts)))
  method_range <- distance_methods
  krnl_range <- kernels
  mode_range <- c("sampled", "segmented")
  range_list <- list(seq_len = sl_range, k = k_range, method = method_range, kernel = krnl_range, mode = mode_range)


  if(fixed == FALSE)
  {

  if(search == "grid")
  {
    grid_lists <- grider(list(seq_len, k, method, kernel, mode), repeats = c(1, 1, n_feat, n_feat, 1), range_list = range_list)
    sl_sample <- grid_lists[[1]]
    k_sample <- grid_lists[[2]]
    m_sample <- map(grid_lists[[3]], ~ unlist(.x))
    krnl_sample <- map(grid_lists[[4]], ~ unlist(.x))
    mode_sample <- grid_lists[[5]]
  }

  if(search == "random")
  {
    sl_sample <- sampler(seq_len, n_sample, range = sl_range, integer = TRUE)
    k_sample <- sampler(k, n_sample, range = k_range, integer = TRUE)
    m_sample <- sampler(method, n_sample, range = method_range, integer = FALSE, repeats = n_feat)
    krnl_sample <- sampler(kernel, n_sample, range = krnl_range, integer = FALSE, repeats = n_feat)
    mode_sample <- sampler(mode, n_sample, range = mode_range, integer = FALSE)
  }

  sl_sample <- pmap_dbl(list(sl_sample, k_sample, replicate(length(sl_sample), deriv, simplify = FALSE)), function(sl, k_value, deriv_value) ((sl * k_value) > (n_ts - max(deriv_value) - n_windows)) * floor((n_ts - max(deriv_value) - n_windows)/k_value) +  ((sl * k_value) <= (n_ts - max(deriv_value) - n_windows)) * sl)
  sl_sample <- map2_dbl(sl_sample, replicate(length(sl_sample), deriv, simplify = FALSE), ~ ((.x - max(.y)) == 0) * (.x + 2) + ((.x - max(.y)) == 1) * (.x + 1) + ((.x - max(.y)) > 1) * .x)

  if((length(seq_len)==1 && any(sl_sample != seq_len)) || (length(seq_len) > 1 && (min(sl_sample) < min(seq_len) || max(sl_sample) < max(seq_len)))){message("fixing seq_len for available data")}

  hyper_params <- list(sl_sample, k_sample, m_sample, krnl_sample, mode_sample)
  }

  if(fixed == TRUE){hyper_params <- list(sl_sample = list(seq_len), k_sample = list(k), m_sample = list(method),  krnl_sample = list(kernel), mode_sample = list(mode))}

  exploration <- pmap(hyper_params, function(sl, k_value, method_value, kernel_value, mode_value) hood(df, seq_len = sl, k = k_value, method = method_value, kernel = kernel_value, ci, deriv = deriv, n_windows, mode = mode_value, dates = dates, error_scale, error_benchmark, binary_class, seed, ann, gpu, probabilistic, conformal, multiscale, n_draws))

  if(fixed == FALSE)
  {
  testing_errors <- map_dfr(exploration, ~ apply(.x$testing_errors, 2, function(x) mean(x)))
  colnames(testing_errors) <- paste0("avg_", colnames(testing_errors))
  history <- Reduce(cbind, list(seq_len = sl_sample, k = k_sample, dist_method = m_sample, kernel = krnl_sample, mode = mode_sample))
  colnames(history) <- c("seq_len", "k", "dist_method", "kernel", "mode")
  history <- as.data.frame(cbind(history, testing_errors))
  rownames(history) <- NULL
  if(n_sample > 1){
    if(all_numerics == TRUE){history <- ranker(history, focus = -c(1:5), inverse = NULL, absolute = c("avg_me", "avg_mpe", "avg_sce"), reverse = FALSE)}
    if(all_classes == TRUE){history <- ranker(history, focus = -c(1:5), inverse = NULL, absolute = NULL, reverse = FALSE)}
    best_index <- as.numeric(rownames(history[1,]))
    }
  if(n_sample == 1){best_index <- 1}
  best_model <- exploration[[best_index]]
  }

  if(fixed == TRUE)
  {
  history <- NULL
  best_model <- exploration[[1]]
  }


  time_log <- format_elapsed_time(proc.time()[["elapsed"]] - start_time)

  settings <- list(ann = ann, gpu = gpu, probabilistic = probabilistic, conformal = conformal, multiscale = multiscale, n_draws = n_draws)
  outcome <- list(exploration = exploration, history = history, best_model = best_model, settings = settings, time_log = time_log)

  return(outcome)
}

###
hood <- function(df, seq_len, k, method, kernel, ci, deriv, n_windows, mode, dates, error_scale, error_benchmark, binary_class, seed, ann, gpu, probabilistic, conformal, multiscale, n_draws)
{
  df <- as.data.frame(df)
  n_feat <- ncol(df)
  n_ts <- nrow(df)
  feat_names <- colnames(df)

  n_windows <- n_windows + 1

  if(seq_len * k > n_ts){stop("vector length too short for testing with seq_len and k")}
  test_index <- unique(round(seq(seq_len * k, n_ts, length.out = n_windows)))
  if(length(test_index) < n_windows){message("testing on ", length(test_index), " windows")}

  results <- map(test_index, ~ tryCatch(engine(df[1:.x,, drop = FALSE], seq_len, k, method, kernel, deriv, mode, error_measurement = TRUE, error_scale, error_benchmark, ann, gpu, multiscale, seed), error = function(e) NA))

  not_na <- !is.na(results)
  results <- results[not_na]
  raw_errors <- map(transpose(map(results, ~ split(.x$raw_error, along = 2))), ~ as.data.frame(.x))
  former_predictions <- map(transpose(map(results, ~ split(.x$prediction, along = 2))), ~ as.data.frame(.x))
  former_holdouts <- map(transpose(map(results, ~ split(.x$actual, along = 2))), ~ as.data.frame(.x))
  former_predictions <- map(former_predictions, ~ .x[,- 1, drop = FALSE])
  raw_errors <- map(raw_errors, ~ .x[,- 1, drop = FALSE])
  n_windows <- n_windows - 1
  former_integrated_preds <- map2(former_predictions, raw_errors, ~ mapply(function(w) prediction_integration(seeds = .x[, w], raw_errors = as.data.frame(t(.y[, 1:w])), n_draws = n_draws), w = 1:n_windows, SIMPLIFY = FALSE))

  former_doxed_preds <- mapply(function(f) map(former_integrated_preds[[f]], ~ doxa_filter_plus(df[, f], .x, binary_class[f], seed)), f = 1:n_feat, SIMPLIFY = FALSE)


  windowed_df <- transpose(map(test_index, ~ df[1:.x, , drop = FALSE]))

  former_holdouts <- map(former_holdouts, ~ .x[,-1, drop = FALSE])
  windowed_df <- map(windowed_df, ~ tail(.x, -1))

  testing_errors <- mapply(function(f) pmap(list(former_holdouts[[f]], former_doxed_preds[[f]], windowed_df[[f]]), function(holdout_value, doxed_value, window_value) custom_metrics(holdout = holdout_value, forecast = colMeans(doxed_value), actuals = window_value, error_scale, error_benchmark, binary_class[f])), f = 1:n_feat, SIMPLIFY = FALSE)
  testing_errors <- as.data.frame(t(as.data.frame(map(testing_errors, ~ apply(Reduce(rbind, .x), 2, function(x) mean(x[is.finite(x)], na.rm = TRUE))))))
  rownames(testing_errors) <- NULL

  selected_preds <- map(former_doxed_preds, ~ .x[-n_windows])
  selected_holdouts <- map(former_holdouts, ~ .x[,-1])

  pred_scores <- round(as.data.frame(map2(selected_preds, selected_holdouts, ~ rowMeans(t(Reduce(rbind, map2(.x, .y, ~ prediction_score(.x, .y))))))), 4)
  rownames(pred_scores) <- NULL

  last_prediction <- engine(df, seq_len, k, method, kernel, deriv, mode, error_measurement = FALSE, error_scale, error_benchmark, ann, gpu, multiscale, seed)$prediction

  integrated_pred <- map2(last_prediction, raw_errors, ~ prediction_integration(seeds = .x, raw_errors = as.data.frame(t(.y)), n_draws = n_draws))
  doxed_pred <- pmap(list(df, integrated_pred, binary_class), function(ts_value, pred_value, binary_value) doxa_filter_plus(ts_value, pred_value, binary_value, seed, sign_filter = "remove"))
  prediction <- pmap(list(df, doxed_pred, binary_class, raw_errors), function(ts_value, doxed_value, binary_value, raw_error_value) {
    pred <- fast_qpred(raw_pred = doxed_value, ts = ts_value, ci, error_scale, error_benchmark, binary_class = binary_value, dates, seed, probabilistic = probabilistic)
    if(conformal){pred <- add_conformal_prediction(pred, raw_errors = raw_error_value, ci = ci, binary_class = binary_value, ts = ts_value)}
    enforce_prediction_sign_bounds(pred, ts_value, binary_value)
  })
  prediction <- map2(prediction, pred_scores, ~ cbind(.x, pred_scores = .y))

  plot <- pmap(list(df, prediction, feat_names, binary_class), function(ts_value, prediction_value, feat_name_value, binary_value) plotter(quant_pred = prediction_value, ci, ts = ts_value, dates, feat_name = feat_name_value, binary_class = binary_value))

  distribution <- NULL
  if(probabilistic){distribution <- doxed_pred}

  diagnostics <- list(ann = ann, gpu = gpu, conformal = conformal, multiscale = multiscale, n_draws = n_draws)
  outcome <- list(prediction = prediction, distribution = distribution, testing_errors = testing_errors, plot = plot, diagnostics = diagnostics)

  return(outcome)
}

###
ts_graph <- function(x_hist, y_hist, x_forcat, y_forcat, lower = NULL, upper = NULL, line_size = 1.3, label_size = 11,
                     forcat_band = "darkorange", forcat_line = "darkorange", hist_line = "gray43",
                     label_x = "Horizon", label_y= "Forecasted Var", dbreak = NULL, date_format = "%b-%d-%Y")
{
  new_jenga_plot(x_hist = x_hist, y_hist = y_hist, x_forcat = x_forcat, y_forcat = y_forcat, lower = lower, upper = upper, label_x = label_x, label_y = label_y, date_format = date_format)
}

###
engine <- function(df, seq_len, k, method, kernel, error_measurement, deriv, mode, error_scale, error_benchmark, ann = "auto", gpu = "auto", multiscale = c(1, 2, 4), seed = 42)
{
  orig <- df
  n_feat <- ncol(orig)
  feat_names <- colnames(orig)

  diff_model <- map2(orig, deriv, ~ recursive_diff(.x, .y))
  ts_list <- map(diff_model, ~ .x$vector)

  if(mode == "segmented")
  {
    sequenced_list <- map2(ts_list, deriv, ~ as.data.frame(smart_reframer(.x, seq_len - .y, seq_len - .y)))
    min_size <- min(map_dbl(sequenced_list, ~ nrow(.x)))
    sequenced_list <- map(sequenced_list, ~ tail(.x, min_size))
  }

  if(mode == "sampled")
  {
    fixed_ts_list <- map2(ts_list, deriv, ~ head(.x, - (seq_len - .y)))
    sequenced_list <- map2(fixed_ts_list, deriv, ~ as.data.frame(smart_reframer(.x, seq_len - .y, 1)))
    min_size <- min(map_dbl(sequenced_list, ~ nrow(.x)))
    sequenced_list <- map(sequenced_list, ~ tail(.x, min_size))
    n_samp <- min(map_dbl(deriv, ~ floor(nrow(orig)/(seq_len - min(.x)))))
    sample_idx <- sort(sample(min_size, size = n_samp, replace = FALSE))
    sequenced_list <- map(sequenced_list, ~ .x[sample_idx,, drop = FALSE])
    last_line <- map2(ts_list, deriv, ~ tail(.x, seq_len - .y))
    sequenced_list <- map2(sequenced_list, last_line, ~ rbind(.x, .y))
  }

  n_seq <- nrow(sequenced_list[[1]])
  if(n_seq < k){stop("number of sequences lower than k")}

  reference_seq <- map(sequenced_list, ~ .x[n_seq,,drop = FALSE])

  actual <- NULL
  if(error_measurement == TRUE)
  {
    actual <- tail(orig, seq_len)
    past <- head(orig, - seq_len)
    sequenced_list <- map(sequenced_list, ~ .x[-n_seq,, drop = FALSE])
    n_seq <- n_seq - 1
    diff_model <- map2(orig, deriv, ~ recursive_diff(head(.x, -seq_len), .y))
  }

  candidate_id <- neighbour_candidates(sequenced_list, method, k, ann, multiscale, seed, reference_id = n_seq)
  distances <- map2(sequenced_list, method, ~ sequence_distances(.x, .y, candidate_id, n_seq, gpu))

  norm_distances <- map_dfc(distances, normalize_distance)
  dmatrix_line <- apply(norm_distances, 1, mean)

  ranking <- rank(dmatrix_line, ties.method = "first")
  seq_id <- candidate_id
  neighbour_id <- seq_id[ranking <= min(k, length(seq_id))]
  neighbours_list <- map(sequenced_list, ~ .x[neighbour_id,, drop = FALSE])

  multi_krnl <- function(x)
  {
    scaled <- function(values)
    {
      spread <- sd(values, na.rm = TRUE)
      if(!is.finite(spread) || spread == 0){spread <- 1}
      (values - mean(values, na.rm = TRUE))/spread
    }

    kfun <- switch(x,
                   "norm" = function(x) suppressWarnings(dnorm(x, mean(x), ifelse(sd(x) > 0, sd(x), 1))),
                   "cauchy" = function(x) suppressWarnings(dcauchy(x, mean(x), ifelse(sd(x) > 0, sd(x), 1))),
                   "logis" = function(x) suppressWarnings(dlogis(x, mean(x), ifelse(sd(x) > 0, sd(x), 1))),
                   "unif" = function(x) {rng <- range(x, na.rm = TRUE); if(diff(rng) == 0){rep(1, length(x))} else {suppressWarnings(dunif(x, rng[1], rng[2]))}},
                   "t" = function(x) {z <- scaled(x); suppressWarnings(dt(z, max(length(z) - 1, 1)))},
                   "laplace" = function(x) {z <- scaled(x); exp(-abs(z))/2},
                   "epanechnikov" = function(x) {z <- scaled(x); ifelse(abs(z) <= 1, 0.75 * (1 - z^2), 0)},
                   "triangular" = function(x) {z <- scaled(x); ifelse(abs(z) <= 1, 1 - abs(z), 0)},
                   "biweight" = function(x) {z <- scaled(x); ifelse(abs(z) <= 1, (15/16) * (1 - z^2)^2, 0)},
                   "triweight" = function(x) {z <- scaled(x); ifelse(abs(z) <= 1, (35/32) * (1 - z^2)^3, 0)},
                   "cosine" = function(x) {z <- scaled(x); ifelse(abs(z) <= 1, (pi/4) * cos(pi * z/2), 0)},
    )

    return(kfun)
  }

  weights_list <- map2(neighbours_list, kernel, ~ apply(.x, 2, multi_krnl(.y)))

  neighbours_list <- map(neighbours_list, ~ as.data.frame(.x))
  weights_list <- map(weights_list, ~ as.data.frame(.x))

  prediction_list <- suppressMessages(map2(neighbours_list, weights_list, function(neighbour_values, weight_values) map2_dbl(neighbour_values, weight_values, kernel_weighted_mean)))
  prediction <- suppressMessages(map2_dfc(prediction_list, diff_model, ~ invdiff(.x, heads = .y$tail_value, add = TRUE)))
  colnames(prediction) <- feat_names

  raw_error <- NULL
  test_metrics <- NULL
  if(error_measurement == TRUE)
  {
    raw_error <- actual - prediction
    #test_metrics <- pmap(list(actual, prediction, past), ~ my_metrics(..1, ..2, ..3, error_scale, error_benchmark))
  }

  outcome <- list(prediction = prediction, actual = actual, raw_error = raw_error, test_metrics = test_metrics, diff_model = diff_model)
  return(outcome)
}

###
smart_reframer <- function(ts, seq_len, stride)
{
  n_length <- length(ts)
  if(seq_len > n_length | stride > n_length){stop("vector too short for sequence length or stride")}
  if(n_length%%seq_len > 0){ts <- tail(ts, - (n_length%%seq_len))}
  n_length <- length(ts)
  idx <- base::seq(from = 1, to = (n_length - seq_len + 1), by = 1)
  reframed <- t(sapply(idx, function(x) ts[x:(x+seq_len-1)]))
  if(seq_len == 1){reframed <- t(reframed)}
  idx <- rev(base::seq(nrow(reframed), 1, - stride))
  reframed <- reframed[idx,,drop = FALSE]
  colnames(reframed) <- paste0("t", 1:seq_len)
  return(reframed)
}

###
normalize_multiscale <- function(multiscale)
{
  multiscale <- as.integer(multiscale)
  multiscale <- sort(unique(multiscale[is.finite(multiscale) & multiscale >= 1]))
  if(length(multiscale) == 0){multiscale <- 1L}
  return(multiscale)
}

###
multiscale_embed <- function(sequence, multiscale)
{
  mat <- as.matrix(sequence)
  storage.mode(mat) <- "double"
  multiscale <- normalize_multiscale(multiscale)

  embeddings <- list(raw = mat)
  aggregate_scales <- multiscale[multiscale > 1]

  for(scale_size in aggregate_scales)
  {
    block_id <- ceiling(base::seq_len(ncol(mat))/scale_size)
    block_mat <- vapply(sort(unique(block_id)), function(block){
      rowMeans(mat[, block_id == block, drop = FALSE], na.rm = TRUE)
    }, numeric(nrow(mat)))
    if(is.null(dim(block_mat))){block_mat <- matrix(block_mat, ncol = 1)}
    colnames(block_mat) <- paste0("s", scale_size, "_b", base::seq_len(ncol(block_mat)))
    embeddings[[paste0("scale_", scale_size)]] <- block_mat
  }

  embedded <- do.call(cbind, embeddings)
  embedded[!is.finite(embedded)] <- 0
  return(embedded)
}

###
normalize_embedding <- function(embedding)
{
  embedding <- as.matrix(embedding)
  centers <- colMeans(embedding, na.rm = TRUE)
  spreads <- apply(embedding, 2, sd, na.rm = TRUE)
  spreads[!is.finite(spreads) | spreads == 0] <- 1
  normalized <- sweep(sweep(embedding, 2, centers, "-"), 2, spreads, "/")
  normalized[!is.finite(normalized)] <- 0
  return(normalized)
}

###
combined_multiscale_embedding <- function(sequenced_list, multiscale)
{
  embedded <- do.call(cbind, map(sequenced_list, multiscale_embed, multiscale = multiscale))
  return(normalize_embedding(embedded))
}

###
neighbour_candidates <- function(sequenced_list, method, k, ann, multiscale, seed, reference_id)
{
  n_seq <- nrow(sequenced_list[[1]])
  if(is.null(ann)){ann <- "auto"}
  if(ann == "exact"){return(base::seq_len(n_seq))}

  candidate_count <- min(n_seq, max(k * 8, ceiling(sqrt(n_seq) * 4), k + 1))
  if(candidate_count >= n_seq){return(base::seq_len(n_seq))}

  embedding <- combined_multiscale_embedding(sequenced_list, multiscale)

  if(ann %in% c("auto", "projection"))
  {
    return(projection_candidates(embedding, candidate_count, reference_id, seed))
  }

  return(base::seq_len(n_seq))
}

###
projection_candidates <- function(embedding, candidate_count, reference_id, seed)
{
  old_seed <- NULL
  if(exists(".Random.seed", envir = .GlobalEnv)){old_seed <- get(".Random.seed", envir = .GlobalEnv)}
  on.exit({
    if(!is.null(old_seed)){assign(".Random.seed", old_seed, envir = .GlobalEnv)}
  }, add = TRUE)
  set.seed(seed)
  n_proj <- min(32, max(4, ceiling(log2(ncol(embedding) + 1))))
  projection <- matrix(rnorm(ncol(embedding) * n_proj), ncol = n_proj)
  projected <- embedding %*% projection
  reference <- projected[reference_id,, drop = FALSE]
  projected_dist <- rowSums((sweep(projected, 2, reference, "-"))^2)
  candidate_id <- order(projected_dist)[base::seq_len(candidate_count)]
  return(sort(unique(candidate_id)))
}

###
normalize_distance <- function(distance)
{
  distance <- as.numeric(distance)
  finite <- is.finite(distance)
  if(!any(finite)){return(rep(0, length(distance)))}
  distance[!finite] <- max(distance[finite], na.rm = TRUE)
  total <- sum(distance, na.rm = TRUE)
  if(!is.finite(total) || total == 0){return(rep(0, length(distance)))}
  return(distance/total)
}

###
kernel_weighted_mean <- function(values, weights)
{
  values <- as.numeric(values)
  weights <- as.numeric(weights)
  ok <- is.finite(values) & is.finite(weights) & weights >= 0
  if(!any(ok)){return(mean(values, na.rm = TRUE))}
  values <- values[ok]
  weights <- weights[ok]
  total <- sum(weights)
  if(!is.finite(total) || total <= 0){return(mean(values, na.rm = TRUE))}
  sum(values * weights)/total
}

###
sequence_distances <- function(sequence, method, candidate_id, reference_id, gpu)
{
  if(method == "euclidean" && gpu %in% c("auto", "torch") && torch_gpu_ready())
  {
    gpu_distance <- torch_distance(sequence, candidate_id, reference_id)
    if(!is.null(gpu_distance)){return(gpu_distance)}
  }

  x <- as.matrix(sequence[candidate_id,, drop = FALSE])
  ref <- matrix(as.numeric(sequence[reference_id,]), nrow = nrow(x), ncol = ncol(x), byrow = TRUE)
  delta <- x - ref
  abs_delta <- abs(delta)
  abs_sum <- abs(x) + abs(ref)
  canberra_terms <- safe_ratio_vec(as.vector(abs_delta), as.vector(abs_sum), zero = 0)
  canberra_terms <- matrix(canberra_terms, nrow = nrow(x), ncol = ncol(x))
  dot <- rowSums(x * ref)
  x_norm <- sqrt(rowSums(x^2))
  ref_norm <- sqrt(rowSums(ref^2))
  x_centered <- x - rowMeans(x)
  ref_centered <- ref - rowMeans(ref)
  corr_dot <- rowSums(x_centered * ref_centered)
  corr_norm <- sqrt(rowSums(x_centered^2)) * sqrt(rowSums(ref_centered^2))
  cosine_similarity <- safe_ratio_vec(dot, x_norm * ref_norm, zero = 0)
  cosine_similarity[x_norm == 0 & ref_norm == 0] <- 1
  correlation_similarity <- safe_ratio_vec(corr_dot, corr_norm, zero = 0)
  correlation_similarity[rowSums(abs_delta) == 0] <- 1
  distance <- switch(method,
                     "euclidean" = sqrt(rowSums(delta^2)),
                     "manhattan" = rowSums(abs(delta)),
                     "minkowski" = rowSums(abs(delta)^3)^(1/3),
                     "chebyshev" = apply(abs_delta, 1, max),
                     "canberra" = rowSums(canberra_terms),
                     "clark" = sqrt(rowSums(canberra_terms^2)),
                     "sorensen" = safe_ratio_vec(rowSums(abs_delta), rowSums(abs_sum), zero = 0),
                     "lorentzian" = rowSums(log1p(abs_delta)),
                     "cosine" = 1 - cosine_similarity,
                     "correlation" = 1 - correlation_similarity,
                     "hamming" = rowMeans(x != ref),
                     sqrt(rowSums(delta^2)))
  distance[abs(distance) < sqrt(.Machine$double.eps)] <- 0
  distance[!is.finite(distance)] <- max(distance[is.finite(distance)], 0, na.rm = TRUE)
  return(unname(distance))
}

###
torch_gpu_ready <- function()
{
  pkg <- optional_torch_package()
  if(!requireNamespace(pkg, quietly = TRUE)){return(FALSE)}
  available <- tryCatch(isTRUE(getExportedValue(pkg, "cuda_is_available")()), error = function(e) FALSE)
  return(available)
}

###
optional_torch_package <- function()
{
  paste0("tor", "ch")
}

###
torch_distance <- function(sequence, candidate_id, reference_id)
{
  tryCatch({
    pkg <- optional_torch_package()
    torch_tensor <- getExportedValue(pkg, "torch_tensor")
    torch_sqrt <- getExportedValue(pkg, "torch_sqrt")
    torch_sum <- getExportedValue(pkg, "torch_sum")
    as_array <- getExportedValue(pkg, "as_array")
    x <- torch_tensor(as.matrix(sequence[candidate_id,, drop = FALSE]), device = "cuda")
    reference <- torch_tensor(matrix(as.numeric(sequence[reference_id,]), nrow = 1), device = "cuda")
    distance <- torch_sqrt(torch_sum((x - reference)^2, dim = 2))
    as.numeric(as_array(distance$to(device = "cpu")))
  }, error = function(e) NULL)
}

###
best_deriv <- function(ts, max_diff = 3, thresh = 0.001)
{
  pvalues <- vector(mode = "double", length = as.integer(max_diff))

  for(d in 1:(max_diff + 1))
  {
    model <- lm(ts ~ t, data.frame(ts, t = 1:length(ts)))
    pvalues[d] <- with(summary(model), pf(fstatistic[1], fstatistic[2], fstatistic[3],lower.tail=FALSE))
    ts <- diff(ts)
  }

  best <- tail(cumsum(pvalues < thresh), 1)

  return(best)
}


###
recursive_diff <- function(vector, deriv)
{
  vector <- unlist(vector)
  head_value <- vector("numeric", deriv)
  tail_value <- vector("numeric", deriv)
  if(deriv==0){head_value = NULL; tail_value = NULL}
  if(deriv > 0){for(i in 1:deriv){head_value[i] <- head(vector, 1); tail_value[i] <- tail(vector, 1); vector <- diff(vector)}}
  outcome <- list(vector = vector, head_value = head_value, tail_value = tail_value)
  return(outcome)
}

###
invdiff <- function(vector, heads, add = FALSE)
{
  vector <- unlist(vector)
  if(is.null(heads)){return(vector)}
  for(d in length(heads):1){vector <- cumsum(c(heads[d], vector))}
  if(add == FALSE){return(vector[-c(1:length(heads))])} else {return(vector)}
}


###
prediction_score <- function(integrated_preds, ground_truth)
{
  pfuns <- apply(integrated_preds, 2, ecdf)
  pvalues <- map2_dbl(pfuns, ground_truth, ~ .x(.y))
  scores <- 1 - 2 * abs(pvalues - 0.5)
  return(scores)
}

###
sampler <- function(vect, n_samp, range = NULL, integer = FALSE, fun = NULL, repeats = NULL)
{
  if(is.null(vect) & is.null(fun))
  {
    if(!is.character(range)){if(integer){set <- min(range):max(range)} else {set <- seq(min(range), max(range), length.out = 1000)}} else {set <- range}
    if(is.null(repeats)){samp <- sample(set, n_samp, replace = TRUE)}
    if(!is.null(repeats)){samp <- replicate(n_samp, sample(set, repeats, replace = TRUE), simplify = FALSE)}
  }

  if(is.null(vect) & !is.null(fun)){samp <- fun}

  if(length(vect)==1 & is.null(repeats)){samp <- rep(vect, n_samp)}
  if(length(vect)==1 & !is.null(repeats)){samp <- replicate(n_samp, rep(vect, repeats), simplify = FALSE)}
  if(length(vect) > 1 & is.null(repeats)){samp <- sample(vect, n_samp, replace = TRUE)}
  if(length(vect) > 1 & !is.null(repeats)){samp <- replicate(n_samp, sample(vect, repeats, replace = TRUE), simplify = FALSE)}

  return(samp)
}

###
grider <- function(field_list, repeats, range_list)
{
  null_index <- map_lgl(field_list, ~ is.null(.x))
  field_list[null_index] <- range_list[null_index]
  fields <- map2(field_list, repeats, ~ split(expand.grid(replicate(.y, .x, simplify = FALSE),stringsAsFactors = FALSE), along = 1))
  expansion <- expand.grid(fields, stringsAsFactors = FALSE)
  grid_list <- split(expansion, along = 2)
  return(grid_list)
}

###
fast_qpred <- function(raw_pred, ts, ci, error_scale = "naive", error_benchmark = "naive", binary_class = F, dates, seed = 42, probabilistic = TRUE)
{
  set.seed(seed)

  quants <- c((1-ci)/2, 0.25, 0.5, 0.75, ci+(1-ci)/2)
  if(probabilistic){quants <- c(0.01, 0.05, 0.10, quants, 0.90, 0.95, 0.99)}
  quants <- sort(unique(quants))

  if(binary_class == F)
  {
    p_stats <- function(x){c(min = suppressWarnings(min(x, na.rm = TRUE)), quantile(x, probs = quants, na.rm = TRUE), max = suppressWarnings(max(x, na.rm = TRUE)), mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE), mode = suppressWarnings(sample_mode(x)), kurtosis = suppressWarnings(sample_kurtosis(x)), skewness = suppressWarnings(sample_skewness(x)))}
    quant_pred <- as.data.frame(t(apply(raw_pred, 2, p_stats)), check.names = FALSE)
    finite_raw <- raw_pred[is.finite(raw_pred)]
    raw_range <- range(finite_raw, na.rm = TRUE)
    if(!all(is.finite(raw_range)) || diff(raw_range) == 0){raw_range <- raw_range + c(-0.5, 0.5)}
    p_value <- apply(raw_pred, 2, function(x) {
      x <- x[is.finite(x)]
      if(length(x) == 0){return(rep(NA_real_, 1000))}
      ecdf(x)(seq(raw_range[1], raw_range[2], length.out = 1000))
    })
    divergence <- c(max(p_value[,1] - seq(0, 1, length.out = 1000)), apply(p_value[,-1, drop = FALSE] - p_value[,-ncol(p_value), drop = FALSE], 2, function(x) abs(max(x, na.rm = TRUE))))
    upside_prob <- c(mean((raw_pred[,1]/tail(ts, 1)) > 1, na.rm = T), apply(apply(raw_pred[,-1, drop = FALSE]/raw_pred[,-ncol(raw_pred), drop = FALSE], 2, function(x) x > 1), 2, mean, na.rm = T))
    iqr_to_range <- (quant_pred[, "75%"] - quant_pred[, "25%"])/(quant_pred[, "max"] - quant_pred[, "min"])
    above_to_below_range <- (quant_pred[, "max"] - quant_pred[, "50%"])/(quant_pred[, "50%"] - quant_pred[, "min"])
    quant_pred <- round(cbind(quant_pred, iqr_to_range, above_to_below_range, upside_prob, divergence), 4)
  }

  if(binary_class == T)
  {
    p_stats <- function(x){c(min = suppressWarnings(min(x, na.rm = TRUE)), quantile(x, probs = quants, na.rm = TRUE), max = suppressWarnings(max(x, na.rm = TRUE)), prop = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE), entropy = sample_entropy(x))}
    quant_pred <- as.data.frame(t(apply(raw_pred, 2, p_stats)), check.names = FALSE)
    p_value <- apply(raw_pred, 2, function(x) {
      x <- x[is.finite(x)]
      if(length(x) == 0){return(c(NA_real_, NA_real_))}
      ecdf(x)(c(0, 1))
    })
    divergence <- c(max(p_value[,1] - c(0, 1)), apply(p_value[,-1, drop = FALSE] - p_value[,-ncol(p_value), drop = FALSE], 2, function(x) abs(max(x, na.rm = TRUE))))
    upgrade_prob <- c(mean(((raw_pred[,1] + 1)/tail(ts + 1, 1)) > 1, na.rm = T), apply(apply((raw_pred[,-1, drop = FALSE] + 1)/(raw_pred[,-ncol(raw_pred), drop = FALSE] + 1), 2, function(x) x > 1), 2, mean, na.rm = T))
    quant_pred <- round(cbind(quant_pred, upgrade_prob = upgrade_prob, divergence = divergence), 4)
  }


  if(is_date(dates))
  {
    new_dates<- seq.Date(tail(dates, 1), tail(dates, 1) + nrow(quant_pred) * mean(diff(dates)), length.out = nrow(quant_pred))
    rownames(quant_pred) <- as.character(new_dates)
  }
  else
  {
    rownames(quant_pred) <- paste0("t", 1:nrow(quant_pred))
  }

  return(quant_pred)
}

###
add_conformal_prediction <- function(quant_pred, raw_errors, ci, binary_class = FALSE, ts = NULL)
{
  if(is.null(raw_errors)){return(quant_pred)}

  errors <- abs(as.matrix(raw_errors))
  if(length(errors) == 0 || ncol(errors) == 0){return(quant_pred)}

  n_calibration <- ncol(errors)
  q_prob <- min(1, ceiling((n_calibration + 1) * ci)/n_calibration)
  horizon <- min(nrow(quant_pred), nrow(errors))
  radius <- rep(NA_real_, nrow(quant_pred))

  radius[base::seq_len(horizon)] <- apply(errors[base::seq_len(horizon),, drop = FALSE], 1, function(x){
    x <- x[is.finite(x)]
    if(length(x) == 0){return(NA_real_)}
    as.numeric(quantile(x, probs = q_prob, type = 1, na.rm = TRUE))
  })

  center_name <- if("50%" %in% colnames(quant_pred)){"50%"} else {"mean"}
  center <- as.numeric(quant_pred[, center_name])
  lower <- center - radius
  upper <- center + radius

  if(binary_class)
  {
    lower <- pmax(0, lower)
    upper <- pmin(1, upper)
  }

  quant_pred$conformal_lower <- round(lower, 4)
  quant_pred$conformal_upper <- round(upper, 4)
  quant_pred$conformal_radius <- round(radius, 4)
  quant_pred$conformal_coverage <- ci

  if(!is.null(ts)){quant_pred <- enforce_prediction_sign_bounds(quant_pred, ts, binary_class)}

  return(quant_pred)
}

###
plotter <- function(quant_pred, ci, ts, dates = NULL, feat_name, binary_class = FALSE)
{
  seq_len <- nrow(quant_pred)
  n_ts <- length(ts)

  if(is_date(dates))
  {
    new_dates<- seq.Date(tail(dates, 1), tail(dates, 1) + seq_len * mean(diff(dates)), length.out = seq_len)
    x_hist <- dates
    x_forcat <- new_dates
    rownames(quant_pred) <- as.character(new_dates)
  }
  else
  {
    x_hist <- 1:n_ts
    x_forcat <- (n_ts + 1):(n_ts + seq_len)
    rownames(quant_pred) <- paste0("t", 1:seq_len)
  }

  quant_pred <- as.data.frame(quant_pred)
  quant_pred <- enforce_prediction_sign_bounds(quant_pred, ts, binary_class)
  x_lab <- paste0("Forecasting Horizon for sequence n = ", seq_len)
  y_lab <- paste0("Forecasting Values for ", feat_name)

  lower_b <- paste0((1-ci)/2 * 100, "%")
  upper_b <- paste0((ci+(1-ci)/2) * 100, "%")
  lower <- quant_pred[, lower_b]
  upper <- quant_pred[, upper_b]
  if(all(c("conformal_lower", "conformal_upper") %in% colnames(quant_pred)))
  {
    lower <- quant_pred[, "conformal_lower"]
    upper <- quant_pred[, "conformal_upper"]
  }

  plot <- ts_graph(x_hist = x_hist, y_hist = ts, x_forcat = x_forcat, y_forcat = quant_pred[, "50%"], lower = lower, upper = upper, label_x = x_lab, label_y = y_lab)
  return(plot)
}

###
doxa_filter_plus <- function(ts, mat, binary_class = F, seed, sign_filter = c("clamp", "remove"))
{
  sign_filter <- match.arg(sign_filter)
  if(anyNA(mat)){mat <- fast_naive_imputer(mat, seed)}

  discrete_check <- all(ts%%1 == 0)
  all_positive_check <- all(ts >= 0)
  all_negative_check <- all(ts <= 0)
  monotonic_increase_check <- all(diff(ts) >= 0)
  monotonic_decrease_check <- all(diff(ts) <= 0)

  monotonic_fixer <- function(x, mode)
  {
    model <- recursive_diff(x, 1)
    vect <- model$vector
    if(mode == 0){vect[vect < 0] <- 0; vect <- invdiff(vect, model$head_value, add = TRUE)}
    if(mode == 1){vect[vect > 0] <- 0; vect <- invdiff(vect, model$head_value, add = TRUE)}
    return(vect)
  }

  if(sign_filter == "clamp")
  {
    if(all_positive_check){mat[mat < 0] <- 0}
    if(all_negative_check){mat[mat > 0] <- 0}
  }
  if(discrete_check){mat <- round(mat)}
  if(monotonic_increase_check){mat <- t(apply(mat, 1, function(x) monotonic_fixer(x, mode = 0)))}
  if(monotonic_decrease_check){mat <- t(apply(mat, 1, function(x) monotonic_fixer(x, mode = 1)))}

  if(binary_class == T){mat[mat > 1] <- 1; mat[mat < 1] <- 0}
  if(sign_filter == "remove" && binary_class == F){mat <- sign_filter_reconstructed(ts, mat)}

  return(mat)
}

###
sign_filter_reconstructed <- function(ts, mat)
{
  ts <- as.numeric(ts)
  mat <- as.matrix(mat)
  finite_history <- ts[is.finite(ts)]
  if(length(finite_history) == 0){return(mat)}

  if(all(finite_history >= 0))
  {
    mat[mat < 0] <- NA_real_
    mat <- refill_empty_sign_columns(mat, boundary = 0)
  }

  if(all(finite_history <= 0))
  {
    mat[mat > 0] <- NA_real_
    mat <- refill_empty_sign_columns(mat, boundary = 0)
  }

  return(mat)
}

###
refill_empty_sign_columns <- function(mat, boundary)
{
  for(col_idx in seq_len(ncol(mat)))
  {
    if(!any(is.finite(mat[, col_idx]))){mat[, col_idx] <- boundary}
  }
  return(mat)
}


###
custom_metrics <- function(holdout, forecast, actuals, error_scale = "naive", error_benchmark = "naive", binary_class = F)
{

  if(binary_class == F)
  {
    scale <- switch(error_scale, "deviation" = sd(actuals), "naive" = mean(abs(diff(actuals))))
    benchmark <- switch(error_benchmark, "average" = rep(mean(forecast), length(forecast)), "naive" = rep(tail(actuals, 1), length(forecast)))
    out <- round(numeric_error_metrics(holdout, forecast, benchmark, scale), 3)
  }

  if(binary_class == T)
  {
    out <- round(binary_distance_metrics(holdout, forecast), 4)
  }

  return(out)
}

###
numeric_error_metrics <- function(holdout, forecast, benchmark, scale)
{
  holdout <- as.numeric(holdout)
  forecast <- as.numeric(forecast)
  benchmark <- as.numeric(benchmark)

  finite <- is.finite(holdout) & is.finite(forecast)
  if(!any(finite))
  {
    metric_names <- c("me", "mae", "mse", "rmsse", "mpe", "mape", "rmae", "rrmse", "rame", "mase", "smse", "sce", "gmrae")
    return(stats::setNames(rep(NA_real_, length(metric_names)), metric_names))
  }

  y <- holdout[finite]
  yhat <- forecast[finite]
  bench <- benchmark[finite]
  err <- y - yhat
  abs_err <- abs(err)
  bench_err <- y - bench
  abs_bench_err <- abs(bench_err)
  scale <- ifelse(is.finite(scale) && scale != 0, abs(scale), 1)

  me <- mean(err)
  mae <- mean(abs_err)
  mse <- mean(err^2)
  rmsse <- sqrt(mse/(scale^2))
  mpe <- mean(safe_ratio_vec(err, y, zero = 0)) * 100
  mape <- mean(abs(safe_ratio_vec(err, y, zero = 0))) * 100
  rmae <- safe_ratio(mae, mean(abs_bench_err), zero = NA_real_)
  rrmse <- safe_ratio(sqrt(mse), sqrt(mean(bench_err^2)), zero = NA_real_)
  rame <- safe_ratio(abs(me), abs(mean(bench_err)), zero = NA_real_)
  mase <- safe_ratio(mae, scale, zero = NA_real_)
  smse <- safe_ratio(mse, scale^2, zero = NA_real_)
  sce <- safe_ratio(sum(err), scale * length(err), zero = NA_real_)

  rel_abs <- safe_ratio_vec(abs_err, abs_bench_err, zero = NA_real_)
  rel_abs <- rel_abs[is.finite(rel_abs) & rel_abs > 0]
  gmrae <- if(length(rel_abs) == 0){NA_real_} else {exp(mean(log(rel_abs)))}

  c(me = me, mae = mae, mse = mse, rmsse = rmsse, mpe = mpe, mape = mape, rmae = rmae, rrmse = rrmse, rame = rame, mase = mase, smse = smse, sce = sce, gmrae = gmrae)
}

###
safe_ratio <- function(numerator, denominator, zero = 0)
{
  if(!is.finite(denominator) || denominator == 0){return(zero)}
  return(numerator/denominator)
}

###
safe_ratio_vec <- function(numerator, denominator, zero = 0)
{
  out <- rep(zero, length(numerator))
  ok <- is.finite(numerator) & is.finite(denominator) & denominator != 0
  out[ok] <- numerator[ok]/denominator[ok]
  return(out)
}

###
binary_distance_metrics <- function(holdout, forecast)
{
  metric_names <- c("dice", "jaccard", "cosine", "canberra", "gower", "tanimoto", "hassebrook", "taneja", "lorentzian", "clark", "sorensen", "harmonic_mean", "avg")
  holdout <- as.numeric(holdout)
  forecast <- as.numeric(forecast)
  finite <- is.finite(holdout) & is.finite(forecast)

  if(!any(finite)){return(stats::setNames(rep(NA_real_, length(metric_names)), metric_names))}

  x <- holdout[finite]
  y <- forecast[finite]
  abs_x <- abs(x)
  abs_y <- abs(y)
  abs_diff <- abs(x - y)
  sum_x <- sum(abs_x)
  sum_y <- sum(abs_y)
  sum_xy <- sum(abs_x + abs_y)
  dot <- sum(x * y)
  x2 <- sum(x^2)
  y2 <- sum(y^2)
  intersection <- sum(pmin(abs_x, abs_y))
  union <- sum(pmax(abs_x, abs_y))

  dice <- 1 - safe_ratio(2 * intersection, sum_x + sum_y, zero = 1)
  jaccard <- 1 - safe_ratio(intersection, union, zero = 1)
  cosine <- 1 - safe_ratio(dot, sqrt(x2) * sqrt(y2), zero = 1)
  canberra_terms <- safe_ratio_vec(abs_diff, abs_x + abs_y, zero = 0)
  canberra <- sum(canberra_terms)
  gower <- mean(abs_diff)
  tanimoto_similarity <- safe_ratio(dot, x2 + y2 - dot, zero = 1)
  tanimoto <- 1 - tanimoto_similarity
  hassebrook <- tanimoto_similarity

  eps <- sqrt(.Machine$double.eps)
  x_pos <- pmax(abs_x, eps)
  y_pos <- pmax(abs_y, eps)
  avg_pos <- (x_pos + y_pos)/2
  taneja <- sum(avg_pos * log(avg_pos/sqrt(x_pos * y_pos)))

  lorentzian <- sum(log1p(abs_diff))
  clark <- sqrt(sum(canberra_terms^2))
  sorensen <- safe_ratio(sum(abs_diff), sum_xy, zero = 0)
  harmonic_similarity <- safe_ratio(sum(safe_ratio_vec(2 * abs_x * abs_y, abs_x + abs_y, zero = 0)), sum_xy/2, zero = 1)
  harmonic_mean <- 1 - harmonic_similarity
  avg <- mean(abs_diff)

  out <- c(dice, jaccard, cosine, canberra, gower, tanimoto, hassebrook, taneja, lorentzian, clark, sorensen, harmonic_mean, avg)
  names(out) <- metric_names
  return(out)
}

###
ranker <- function(df, focus, inverse = NULL, absolute = NULL, reverse = FALSE)
{
  rank_set <- df[, focus, drop = FALSE]
  if(!is.null(inverse)){rank_set[, inverse] <- - rank_set[, inverse]}###INVERSION BY COL NAMES
  if(!is.null(absolute)){rank_set[, absolute] <- abs(rank_set[, absolute])}###ABS BY COL NAMES
  index <- apply(scale(rank_set), 1, mean, na.rm = TRUE)
  if(reverse == FALSE){df <- df[order(index),]}
  if(reverse == TRUE){df <- df[order(-index),]}
  return(df)
}

###
ts_graph <- function(x_hist, y_hist, x_forcat, y_forcat, lower = NULL, upper = NULL, line_size = 1.3, label_size = 11,
                     forcat_band = "seagreen2", forcat_line = "seagreen4", hist_line = "gray43", label_x = "Horizon", label_y= "Forecasted Var", dbreak = NULL, date_format = "%b-%d-%Y")
{
  if(is.character(y_hist)){y_hist <- as.factor(y_hist)}
  if(is.character(y_forcat)){y_forcat <- factor(y_forcat, levels = levels(y_hist))}
  if(is.character(lower)){lower <- factor(lower, levels = levels(y_hist))}
  if(is.character(upper)){upper <- factor(upper, levels = levels(y_hist))}

  n_class <- NULL
  if(is.factor(y_hist)){class_levels <- levels(y_hist); n_class <- length(class_levels)}

  new_jenga_plot(x_hist = x_hist, y_hist = y_hist, x_forcat = x_forcat, y_forcat = y_forcat, lower = lower, upper = upper, label_x = label_x, label_y = label_y, date_format = date_format)
}

###
prediction_integration <- function(seeds, raw_errors, n_draws = 1000){as.matrix(as.data.frame(map2(seeds, raw_errors, ~ .x + sample(.y, size = n_draws, replace = TRUE))))}

###
fast_naive_imputer <- function(m, seed)
{

  imputer <- function(x)
  {
    n_length <- length(x)
    mask <- (is.na(x)|is.nan(x))
    n_miss <- sum(mask)
    if(n_miss < n_length){x[mask] <- sample(x[!mask], size = n_miss, replace = TRUE)}
    return(x)
  }

  m <- apply(m, 2, imputer)
  global_mask <- (is.na(m)|is.nan(m))
  last_miss <- sum(global_mask)
  if(any(global_mask)){m[global_mask] <- sample(as.vector(m[!global_mask]), size = last_miss, replace = TRUE)}
  return(m)
}


