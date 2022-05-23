#' jenga
#'
#' @param ts A data frame with time features on columns
#' @param seq_len Positive integer. Time-step number of the projected sequence
#' @param k Positive integer. Number of neighbors to consider when applying kernel average. Min number is 3. Default: NULL (automatic selection).
#' @param method Positive integer. Distance method for calculating neighbors. Possibile options are: "euclidean", "manhattan", "canberra1", "minimum", "maximum", "minkowski", "bhattacharyya", "kullback_leibler", "jensen_shannon". Default: NULL (automatic selection).
#' @param kernel String. Distribution used to calculate kernel densities. Possible options are: "norm", "cauchy", "logis", "unif", "t". Default: NULL (automatic selection).
#' @param ci Confidence interval. Default: 0.8
#' @param n_windows Positive integer. Number of validation tests to measure/sample error. Default: 10.
#' @param mode String. Sequencing method: deterministic ("segmented"), or non-deterministic ("sampled"). Default: NULL (automatic selection).
#' @param n_sample Positive integer. Number of samples for grid or random search. Default: 30.
#' @param search String. Two option available: "grid", "random". Default: "random".
#' @param dates Date. Vector with dates for time features.
#' @param error_scale String. Scale for the scaled error metrics. Two options: "naive" (average of naive one-step absolute error for the historical series) or "deviation" (standard error of the historical series). Default: "naive".
#' @param error_benchmark String. Benchmark for the relative error metrics. Two options: "naive" (sequential extension of last value) or "average" (mean value of true sequence). Default: "naive".
#' @param seed Positive integer. Random seed. Default: 42.
#'
#' @author Giancarlo Vercellino \email{giancarlo.vercellino@gmail.com}
#'
#' @return This function returns a list including:
#' \itemize{
#' \item exploration: list of all not-null models, complete with predictions, test metrics, prediction stats and plot
#' \item history: a table with the sampled models, hyper-parameters, validation errors, weighted average rank
#' \item best_model: results for the best model in term of weighted average rank, including:
#' \itemize{
#' \item predictions: min, max, q25, q50, q75, quantiles at selected ci, mean, sd, mode, skewness, kurtosis, IQR to range, risk ratio, upside probability and divergence for each point fo predicted sequences
#' \item testing_errors: training and testing errors for one-step and sequence for each ts feature (me, mae, mse, rmsse, mpe, mape, rmae, rrmse, rame, mase, smse, sce, gmrae)
#' \item pred_scores: a measure of prediction interval fit for each point in predicted sequence (value range from 0, out of boundaries, to 1, close to the median)
#' }
#' \item time_log
#' }
#'
#' @export
#'
#' @import purrr
#' @import tictoc
#' @importFrom readr parse_number
#' @importFrom lubridate seconds_to_period is.Date as.duration
#' @importFrom stats weighted.mean lm pnorm rnorm
#' @importFrom imputeTS na_kalman
#' @importFrom narray split
#' @importFrom modeest mlv1
#' @importFrom moments skewness kurtosis
#' @importFrom stats quantile sd ecdf
#' @importFrom scales number
#' @importFrom lubridate is.Date
#' @importFrom narray split
#' @importFrom utils head tail
#' @importFrom modeest mlv1
#' @importFrom moments skewness kurtosis
#' @importFrom stats dcauchy density dexp dlnorm dlogis dnorm dt dunif median sd
#' @import dplyr
#' @import abind
#' @import philentropy
#' @import abind
#' @import ggplot2
#' @import greybox
#' @importFrom Rfast Dist
#' @importFrom purrr map map_dbl pmap_dbl map2_dbl pmap map_dfr map_depth transpose map2



#'@examples
#'jenga(covid_in_europe[, c(2, 3)], n_sample = 1)
#'jenga(covid_in_europe[, c(4, 5)], n_sample = 1)
#'
#'

jenga <- function(ts, seq_len = NULL, k = NULL, method = NULL, kernel = NULL, ci = 0.8, n_windows = 10, mode = NULL, n_sample = 30, search = "random", dates = NULL, error_scale = "naive", error_benchmark = "naive", seed = 42)
{
  set.seed(seed)

  tic.clearlog()
  tic("time")

  if(!is.data.frame(ts)){stop("time features must be in dataframe format")}
  if(anyNA(ts)){ts <- as.data.frame(na_kalman(ts))}

  n_ts <- nrow(ts)
  n_feat <- ncol(ts)
  deriv <- map_dbl(ts, ~ best_deriv(.x))
  fixed <- FALSE
  if(all(c(length(seq_len)==1, length(k)==1, length(method)==n_feat, length(kernel)==n_feat, length(mode)==1))){fixed <- TRUE}

  if(!is.null(k) && any(k < 3)){k[k < 3] <- 3; message("setting minimum k value to 3")}
  if(!is.null(seq_len) && any(seq_len < 3)){seq_len[seq_len < 3] <- 3; message("setting minimum seq_len value to 3")}

  ###SAMPLERS
  sl_range <- c(3:floor(sqrt(n_ts)))
  k_range <- c(3:floor(sqrt(n_ts)))
  method_range <- c("euclidean", "manhattan", "canberra1", "minimum", "maximum", "minkowski", "bhattacharyya", "kullback_leibler", "jensen_shannon")
  krnl_range <- c("norm", "cauchy", "logis", "unif", "t")
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
    k_sample <- sampler(seq_len, n_sample, range = k_range, integer = TRUE)
    m_sample <- sampler(method, n_sample, range = method_range, integer = FALSE, repeats = n_feat)
    krnl_sample <- sampler(kernel, n_sample, range = krnl_range, integer = FALSE, repeats = n_feat)
    mode_sample <- sampler(mode, n_sample, range = mode_range, integer = FALSE)
  }

  sl_sample <- pmap_dbl(list(sl_sample, k_sample, replicate(length(sl_sample), deriv, simplify = FALSE)), ~ ((..1 * ..2) > (n_ts - max(..3) - n_windows)) * floor((n_ts - max(..3) - n_windows)/..2) +  ((..1 * ..2) <= (n_ts - max(..3) - n_windows)) * ..1)
  sl_sample <- map2_dbl(sl_sample, replicate(length(sl_sample), deriv, simplify = FALSE), ~ ((.x - max(.y)) == 0) * (.x + 2) + ((.x - max(.y)) == 1) * (.x + 1) + ((.x - max(.y)) > 1) * .x)

  if(length(seq_len)==1 && any(sl_sample != seq_len)){message("fixing seq_len for available data")}

  hyper_params <- list(sl_sample, k_sample, m_sample, krnl_sample, mode_sample)
  }

  if(fixed == TRUE){hyper_params <- list(sl_sample = list(seq_len), k_sample = list(k), m_sample = list(method),  krnl_sample = list(kernel), mode_sample = list(mode))}

  exploration <- suppressMessages(pmap(hyper_params, ~ tryCatch(suppressWarnings(hood(ts, seq_len = ..1, k = ..2, method = ..3, kernel = ..4, ci, deriv = deriv, n_windows, mode = ..5, dates = dates, error_scale, error_benchmark, seed)), error = function(e) NA)))
  not_na <- !is.na(exploration)
  exploration <- exploration[not_na]


  if(fixed == FALSE)
  {
  testing_errors <- map_dfr(exploration, ~ apply(.x$testing_errors, 2, function(x) mean(x)))
  colnames(testing_errors) <- paste0("avg_", colnames(testing_errors))
  history <- Reduce(cbind, list(seq_len = sl_sample[not_na], k = k_sample[not_na], dist_method = m_sample[not_na], kernel = krnl_sample[not_na], mode = mode_sample[not_na]))
  colnames(history) <- c("seq_len", "k", "dist_method", "kernel", "mode")
  history <- as.data.frame(cbind(history, testing_errors))
  rownames(history) <- NULL
  if(n_sample > 1){history <- history[order(unlist(history$avg_mae), unlist(history$avg_mse)),]
  best_index <- as.numeric(rownames(history[1,]))}
  if(n_sample == 1){best_index <- 1}
  best_model <- exploration[[best_index]]
  }

  if(fixed == TRUE)
  {
  history <- NULL
  best_model <- exploration[[1]]
  }


  toc(log = TRUE)
  time_log <- seconds_to_period(round(parse_number(unlist(tic.log())), 0))

  outcome <- list(exploration = exploration, history = history, best_model = best_model, time_log = time_log)

  return(outcome)
}

###
hood <- function(ts, seq_len, k, method, kernel, ci, deriv, n_windows, mode, dates, error_scale, error_benchmark, seed)
{
  ts <- as.data.frame(ts)
  n_feat <- ncol(ts)
  n_ts <- nrow(ts)
  feat_names <- colnames(ts)

  if(seq_len * k > n_ts){stop("vector length too short for testing with seq_len and k")}
  test_index <- unique(round(seq(seq_len * k, n_ts, length.out = n_windows)))
  if(length(test_index) < n_windows){message("testing on ", length(test_index), " windows")}

  results <- map(test_index, ~ tryCatch(engine(ts[1:.x,, drop = FALSE], seq_len, k, method, kernel, deriv, mode, error_measurement = TRUE, error_scale, error_benchmark), error = function(e) NA))

  not_na <- !is.na(results)
  results <- results[not_na]
  raw_errors <- map(transpose(map(results, ~ split(.x$raw_error, along = 2))), ~ as.data.frame(.x))
  former_predictions <- map(transpose(map(results, ~ split(.x$prediction, along = 2))), ~ as.data.frame(.x))
  former_holdouts <- map(transpose(map(results, ~ split(.x$actual, along = 2))), ~ as.data.frame(.x))

  quantile_processing <- pmap(list(former_predictions, raw_errors, former_holdouts, ts), ~
  mapply(function(i) quantile_prediction(raw_pred = NULL, seed_pred = ..1[,i], raw_error = as.data.frame(t(..2[, 1:i])), holdout = ..3[, i], tail_value = NULL, ts = ..4, ci, error_scale, error_benchmark, seed), i = 1: n_windows, SIMPLIFY = FALSE))

  testing_errors <- as.data.frame(t(as.data.frame(map(map_depth(quantile_processing, 2, ~ .x$testing_error), ~ apply(Reduce(rbind, .x), 2, function(x) mean(x[is.finite(x)], na.rm = TRUE))))))
  selected_preds <- map(map_depth(quantile_processing, 2, ~ .x$raw_pred), ~ .x[-n_windows])
  selected_holdouts <- map(former_holdouts, ~ .x[,-1])

  pred_scores <- round(as.data.frame(map2(selected_preds, selected_holdouts, ~ rowMeans(t(Reduce(rbind, map2(.x, .y, ~ prediction_score(.x, .y))))))), 4)

  last_prediction <- engine(ts, seq_len, k, method, kernel, deriv, mode, error_measurement = FALSE, error_scale, error_benchmark)$prediction
  prediction <- map(pmap(list(last_prediction, raw_errors, ts), ~ quantile_prediction(raw_pred = NULL, seed_pred = ..1, raw_error = as.data.frame(t(..2)), holdout = NULL, tail_value = NULL, ts = ..3, ci, error_scale, error_benchmark, seed)), ~ .x$quant_pred)
  prediction <- map(prediction, ~ {rownames(.x) <- paste0("t", 1:seq_len); return(.x)})

  x_hist <- 1:n_ts
  x_forcat <- (n_ts + 1):(n_ts + seq_len)
  rownames(pred_scores) <- paste0("t", 1:seq_len)

  if(is.Date(dates) & length(dates) == n_ts)
  {
    pred_dates <- tail(dates, seq_len) + seq_len
    prediction <- map(prediction, ~ {rownames(.x) <- pred_dates; return(.x)})
    x_hist <- dates
    x_forcat <- pred_dates
    rownames(pred_scores) <- pred_dates
  }

  plot <- pmap(list(prediction, ts, feat_names), ~ ts_graph(x_hist = x_hist, y_hist = ..2, x_forcat = x_forcat, y_forcat = ..1[, 3], lower = ..1[, 1], upper = ..1[, 5], label_y = ..3))

  outcome <- list(prediction = prediction, testing_errors = testing_errors, pred_scores = pred_scores, plot = plot)

  return(outcome)
}

###
ts_graph <- function(x_hist, y_hist, x_forcat, y_forcat, lower = NULL, upper = NULL, line_size = 1.3, label_size = 11,
                     forcat_band = "darkorange", forcat_line = "darkorange", hist_line = "gray43",
                     label_x = "Horizon", label_y= "Forecasted Var", dbreak = NULL, date_format = "%b-%d-%Y")
{
  all_data <- data.frame(x_all = c(x_hist, x_forcat), y_all = c(y_hist, y_forcat))
  forcat_data <- data.frame(x_forcat = x_forcat, y_forcat = y_forcat)

  if(!is.null(lower) & !is.null(upper)){forcat_data$lower <- lower; forcat_data$upper <- upper}

  plot <- ggplot()+geom_line(data = all_data, aes_string(x = "x_all", y = "y_all"), color = hist_line, size = line_size)
  if(!is.null(lower) & !is.null(upper)){plot <- plot + geom_ribbon(data = forcat_data, aes(x = x_forcat, ymin = lower, ymax = upper), alpha = 0.3, fill = forcat_band)}
  plot <- plot + geom_line(data = forcat_data, aes(x = x_forcat, y = y_forcat), color = forcat_line, size = line_size)
  if(!is.null(dbreak)){plot <- plot + scale_x_date(name = paste0("\n", label_x), date_breaks = dbreak, date_labels = date_format)}
  if(is.null(dbreak)){plot <- plot + xlab(label_x)}
  plot <- plot + scale_y_continuous(name = paste0(label_y, "\n"), labels = number)
  plot <- plot + ylab(label_y)  + theme_bw()
  plot <- plot + theme(axis.text=element_text(size=label_size), axis.title=element_text(size=label_size + 2))

  return(plot)
}

###
engine <- function(ts, seq_len, k, method, kernel, error_measurement, deriv, mode, error_scale, error_benchmark)
{
  orig <- ts
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

  distances <- map2(sequenced_list, method, ~ mapply(function(x) Dist(rbind(.x[n_seq,], .x[x,]), method = .y, p = 3), x = 1:n_seq))

  norm_distances <- map_dfc(distances, ~ .x/sum(.x))
  dmatrix_line <- apply(norm_distances, 1, mean)

  ranking <- rank(dmatrix_line)
  seq_id <- 1:n_seq
  neighbour_id <- seq_id[ranking <= k]
  neighbours_list <- map(sequenced_list, ~ .x[neighbour_id,, drop = FALSE])

  multi_krnl <- function(x)
  {
    kfun <- switch(x,
                   "norm" = function(x) suppressWarnings(dnorm(x, mean(x), sd(x))),
                   "cauchy" = function(x) suppressWarnings(dcauchy(x, mean(x), sd(x))),
                   "logis" = function(x) suppressWarnings(dlogis(x, mean(x), sd(x))),
                   "unif" = function(x) suppressWarnings(dunif(x, min(x), max(x))),
                   "t" = function(x) suppressWarnings(dt(x, length(x) - 1)),
    )

    return(kfun)
  }

  weights_list <- map2(neighbours_list, kernel, ~ apply(.x, 2, multi_krnl(.y)))

  neighbours_list <- map(neighbours_list, ~ as.data.frame(.x))
  weights_list <- map(weights_list, ~ as.data.frame(.x))

  prediction_list <- suppressMessages(map2(neighbours_list, weights_list, ~ map2_dbl(.x, .y, ~  sum(..1 * ..2)/sum(..2))))
  prediction <- suppressMessages(map2_dfc(prediction_list, diff_model, ~ invdiff(.x, heads = .y$tail_value, add = TRUE)))
  colnames(prediction) <- feat_names

  raw_error <- NULL
  test_metrics <- NULL
  if(error_measurement == TRUE)
  {
    raw_error <- actual - prediction
    test_metrics <- pmap(list(actual, prediction, past), ~ my_metrics(..1, ..2, ..3, error_scale, error_benchmark))
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
quantile_prediction <- function(raw_pred = NULL, seed_pred = NULL, raw_error = NULL, holdout = NULL, tail_value = NULL, ts, ci, error_scale = "naive", error_benchmark = "naive", seed = 42)
{
  set.seed(seed)

  not_null<- c(!is.null(seed_pred), !is.null(raw_error))
  if(all(not_null))
  {
    pred_integration <- function(pred_seed, error_raw){as.data.frame(map2(pred_seed, error_raw, ~ .x + sample(.y, size = 1000, replace = TRUE)))}
    raw_pred <- t(pred_integration(seed_pred, raw_error))
    if(!is.null(tail_value)){raw_pred <- apply(raw_pred, 2, function(x) invdiff(x, heads = tail_value, add = FALSE))}
  }

  if(!is.null(raw_pred))
  {
    raw_pred <- doxa_filter(ts, t(raw_pred))
    quants <- sort(unique(c((1-ci)/2, 0.25, 0.5, 0.75, ci+(1-ci)/2)))
    p_stats <- function(x){stats <- c(min = suppressWarnings(min(x, na.rm = TRUE)), quantile(x, probs = quants, na.rm = TRUE), max = suppressWarnings(max(x, na.rm = TRUE)), mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE), mode = tryCatch(suppressWarnings(mlv1(x[is.finite(x)], method = "shorth")), error = function(e) NA), kurtosis = tryCatch(suppressWarnings(kurtosis(x[is.finite(x)], na.rm = TRUE)), error = function(e) NA), skewness = tryCatch(suppressWarnings(skewness(x[is.finite(x)], na.rm = TRUE)), error = function(e) NA)); return(stats)}
    quant_pred <- as.data.frame(t(as.data.frame(apply(raw_pred, 2, p_stats))))
    iqr_to_range <- tryCatch((quant_pred[, "75%"] - quant_pred[, "25%"])/(quant_pred[, "max"] - quant_pred[, "min"]), error = function(e) NA)
    median_range_ratio <- tryCatch((quant_pred[, "max"] - quant_pred[, "50%"])/(quant_pred[, "50%"] - quant_pred[, "min"]), error = function(e) NA)
    growths <- mapply(function(m, s) rnorm(1000, m, s), m = quant_pred[, "mean"], s = quant_pred[, "sd"])
    upside_prob <- tryCatch(c(NA, colMeans(apply(growths[,-1, drop = FALSE]/growths[,-ncol(growths), drop = FALSE], 2, function(x) x > 1))), error = function(e) NA)
    pvalue <- mapply(function(m, s) pnorm(seq(min(quant_pred[, "min"]), max(quant_pred[, "max"]), length.out = 1000), m, s), m = quant_pred[, "mean"], s = quant_pred[, "sd"])
    divergence <- tryCatch(c(NA, apply(pvalue[,-1, drop = FALSE] - pvalue[,-ncol(pvalue), drop = FALSE], 2, function(x) abs(max(x, na.rm = TRUE)))), error = function(e) NA)
    quant_pred <- round(cbind(quant_pred, iqr_to_range, median_range_ratio, upside_prob, divergence), 4)
  }

  testing_error <- NULL
  if(!is.null(holdout))
  {
    mean_pred <- colMeans(raw_pred)
    testing_error <- my_metrics(holdout, mean_pred, ts, error_scale, error_benchmark)
  }

  outcome <- list(raw_pred = raw_pred, quant_pred = quant_pred, testing_error = testing_error)
  return(outcome)
}

###
doxa_filter <- function(ts, mat, n_class = NULL)
{
  dim_orig <- dim(mat)
  discrete_check <- all(ts%%1 == 0)
  all_positive_check <- all(ts >= 0)
  all_negative_check <- all(ts <= 0)
  monotonic_increase_check <- all(diff(ts) >= 0)
  monotonic_decrease_check <- all(diff(ts) <= 0)
  class_check <- FALSE
  if(is.integer(n_class)){class_check <- length(unique(ts)) <= n_class}

  monotonic_fixer <- function(x, mode)
  {
    model <- recursive_diff(x, 1)
    vect <- model$vector
    if(mode == 0){vect[vect < 0] <- 0; vect <- invdiff(vect, model$head_value, add = TRUE)}
    if(mode == 1){vect[vect > 0] <- 0; vect <- invdiff(vect, model$head_value, add = TRUE)}
    return(vect)
  }

  if(discrete_check){mat <- floor(mat)}
  if(monotonic_increase_check){mat <- t(apply(mat, 1, function(x) monotonic_fixer(x, mode = 0)))}
  if(monotonic_decrease_check){mat <- t(apply(mat, 1, function(x) monotonic_fixer(x, mode = 1)))}
  if(class_check){mat[!(mat %in% unique(ts))] <- ((mat[!(mat %in% unique(ts))] > max(unique(ts))) * max(unique(ts))) + ((mat[!(mat %in% unique(ts))] < min(unique(ts))) * min(unique(ts)))}
  if(all_positive_check){mat[mat < 0] <- 0}
  if(all_negative_check){mat[mat > 0] <- 0}

  dim(mat) <- dim_orig
  return(mat)
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
my_metrics <- function(holdout, forecast, actuals, error_scale = "naive", error_benchmark = "naive")
{
  scale <- switch(error_scale, "deviation" = sd(actuals), "naive" = mean(abs(diff(actuals))))
  benchmark <- switch(error_benchmark, "average" = rep(mean(forecast), length(forecast)), "naive" = rep(tail(actuals, 1), length(forecast)))

  me <- round(ME(holdout, forecast, na.rm = TRUE), 4)
  mae <- round(MAE(holdout, forecast, na.rm = TRUE), 4)
  mse <- round(MSE(holdout, forecast, na.rm = TRUE), 4)
  rmsse <- round(RMSSE(holdout, forecast, scale, na.rm = TRUE), 4)
  mre <- round(MRE(holdout, forecast, na.rm = TRUE), 4)
  mpe <- round(MPE(holdout, forecast, na.rm = TRUE), 4)
  mape <- round(MAPE(holdout, forecast, na.rm = TRUE), 4)
  rmae <- round(rMAE(holdout, forecast, benchmark, na.rm = TRUE), 4)
  rrmse <- round(rRMSE(holdout, forecast, benchmark, na.rm = TRUE), 4)
  rame <- round(rAME(holdout, forecast, benchmark, na.rm = TRUE), 4)
  mase <- round(MASE(holdout, forecast, scale, na.rm = TRUE), 4)
  smse <- round(sMSE(holdout, forecast, scale, na.rm = TRUE), 4)
  sce <- round(sCE(holdout, forecast, scale, na.rm = TRUE), 4)
  gmrae <- round(GMRAE(holdout, forecast, benchmark, na.rm = TRUE), 4)

  out <- c(me = me, mae = mae, mse = mse, rmsse = rmsse, mpe = mpe, mape = mape, rmae = rmae, rrmse = rrmse, rame = rame, mase = mase, smse = smse, sce = sce, gmrae = gmrae)
  return(out)
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
    samp <- sample(set, n_samp, replace = TRUE)
  }

  if(is.null(vect) & !is.null(fun)){samp <- fun}

  if(length(vect)==1){samp <- rep(vect, n_samp)}
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


