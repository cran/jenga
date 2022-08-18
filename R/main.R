#' jenga
#'
#' @param df A data frame with time features on columns (numerical or categorical features, but not both).
#' @param seq_len Positive integer. Time-step number of the projected sequence
#' @param smoother Logical. Perform optimal smoothing using standard loess (only for numerical features). Default: FALSE
#' @param k Positive integer. Number of neighbors to consider when applying kernel average. Min number is 3. Default: NULL (automatic selection).
#' @param method Positive integer. Distance method for calculating neighbors. Possibile options are: "euclidean", "manhattan", "minkowski". Default: NULL (automatic selection).
#' @param kernel String. Distribution used to calculate kernel densities. Possible options are: "norm", "cauchy", "unif", "t". Default: NULL (automatic selection).
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
#' \item exploration: list of all models, complete with predictions, test metrics, prediction stats and plot
#' \item history: a table with the sampled models, hyper-parameters, validation errors
#' \item best_model: results for the best model, including:
#' \itemize{
#' \item predictions: min, max, q25, q50, q75, quantiles at selected ci, and different statics for numerical and categorical variables
#' \item testing_errors: training and testing errors for one-step and sequence for each ts feature (different measures for numerical and categorical variables)
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
#' @importFrom entropy entropy
#' @importFrom fANCOVA loess.as
#' @import fastDummies


#'@examples
#'jenga(covid_in_europe[, c(2, 3)], n_sample = 1)
#'jenga(covid_in_europe[, c(4, 5)], n_sample = 1)
#'
#'

jenga <- function(df, seq_len = NULL, smoother = FALSE, k = NULL, method = NULL, kernel = NULL, ci = 0.8, n_windows = 10, mode = NULL, n_sample = 30, search = "random", dates = NULL, error_scale = "naive", error_benchmark = "naive", seed = 42)
{
  set.seed(seed)

  tic.clearlog()
  tic("time")


  if(!is.data.frame(df)){stop("time features must be in dataframe format")}

  n_ts <- nrow(df)

  class_index <- any(map_lgl(df, ~ is.factor(.x) | is.character(.x)))
  all_classes <- all(class_index)
  numeric_index <- map_lgl(df, ~ is.integer(.x) | is.numeric(.x))
  all_numerics <- all(numeric_index)
  if(!(all_classes | all_numerics)){stop("only all numerics or all classes, not both")}

  if(all_classes){df <- dummy_cols(df, select_columns = NULL, remove_first_dummy = FALSE, remove_most_frequent_dummy = TRUE, ignore_na = FALSE, split = NULL, remove_selected_columns = TRUE); binary_class <- rep(T, ncol(df))}
  if(all_numerics){binary_class <- rep(FALSE, ncol(df))}

  if(anyNA(df) & all_numerics){df <- as.data.frame(na_kalman(df)); message("kalman imputation on time features\n")}
  if(anyNA(df) & all_classes){df <- floor(as.data.frame(na_kalman(df))); message("kalman imputation on time features\n")}
  if(smoother == TRUE & all_numerics){df <- as.data.frame(purrr::map(df, ~ suppressWarnings(loess.as(x=1:n_length, y=.x)$fitted))); message("performing optimal smoothing\n")}

  n_feat <- ncol(df)
  feat_names <- colnames(df)

  deriv <- map_dbl(df, ~ best_deriv(.x))
  fixed <- FALSE
  if(all(c(length(seq_len)==1, length(k)==1, length(method)==n_feat, length(kernel)==n_feat, length(mode)==1))){fixed <- TRUE}

  if(!is.null(k) && any(k < 3)){k[k < 3] <- 3; message("setting minimum k value to 3")}
  if(!is.null(seq_len) && any(seq_len < 3)){seq_len[seq_len < 3] <- 3; message("setting minimum seq_len value to 3")}

  ###SAMPLERS
  sl_range <- c(3:floor(sqrt(n_ts)))
  k_range <- c(3:floor(sqrt(n_ts)))
  method_range <- c("euclidean", "manhattan", "minkowski")###OUT:  "canberra1", "minimum", "maximum", "bhattacharyya", "kullback_leibler", "jensen_shannon"
  krnl_range <- c("norm", "cauchy", "unif", "t")###OUT: LOGIS
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

  sl_sample <- pmap_dbl(list(sl_sample, k_sample, replicate(length(sl_sample), deriv, simplify = FALSE)), ~ ((..1 * ..2) > (n_ts - max(..3) - n_windows)) * floor((n_ts - max(..3) - n_windows)/..2) +  ((..1 * ..2) <= (n_ts - max(..3) - n_windows)) * ..1)
  sl_sample <- map2_dbl(sl_sample, replicate(length(sl_sample), deriv, simplify = FALSE), ~ ((.x - max(.y)) == 0) * (.x + 2) + ((.x - max(.y)) == 1) * (.x + 1) + ((.x - max(.y)) > 1) * .x)

  if((length(seq_len)==1 && any(sl_sample != seq_len)) || (length(seq_len) > 1 && (min(sl_sample) < min(seq_len) || max(sl_sample) < max(seq_len)))){message("fixing seq_len for available data")}

  hyper_params <- list(sl_sample, k_sample, m_sample, krnl_sample, mode_sample)
  }

  if(fixed == TRUE){hyper_params <- list(sl_sample = list(seq_len), k_sample = list(k), m_sample = list(method),  krnl_sample = list(kernel), mode_sample = list(mode))}

  exploration <- pmap(hyper_params, ~ hood(df, seq_len = ..1, k = ..2, method = ..3, kernel = ..4, ci, deriv = deriv, n_windows, mode = ..5, dates = dates, error_scale, error_benchmark, binary_class, seed))

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


  toc(log = TRUE)
  time_log <- seconds_to_period(round(parse_number(unlist(tic.log())), 0))

  outcome <- list(exploration = exploration, history = history, best_model = best_model, time_log = time_log)

  return(outcome)
}

###
hood <- function(df, seq_len, k, method, kernel, ci, deriv, n_windows, mode, dates, error_scale, error_benchmark, binary_class, seed)
{
  df <- as.data.frame(df)
  n_feat <- ncol(df)
  n_ts <- nrow(df)
  feat_names <- colnames(df)

  n_windows <- n_windows + 1

  if(seq_len * k > n_ts){stop("vector length too short for testing with seq_len and k")}
  test_index <- unique(round(seq(seq_len * k, n_ts, length.out = n_windows)))
  if(length(test_index) < n_windows){message("testing on ", length(test_index), " windows")}

  results <- map(test_index, ~ tryCatch(engine(df[1:.x,, drop = FALSE], seq_len, k, method, kernel, deriv, mode, error_measurement = TRUE, error_scale, error_benchmark), error = function(e) NA))

  not_na <- !is.na(results)
  results <- results[not_na]
  raw_errors <- map(transpose(map(results, ~ split(.x$raw_error, along = 2))), ~ as.data.frame(.x))
  former_predictions <- map(transpose(map(results, ~ split(.x$prediction, along = 2))), ~ as.data.frame(.x))
  former_holdouts <- map(transpose(map(results, ~ split(.x$actual, along = 2))), ~ as.data.frame(.x))
  former_predictions <- map(former_predictions, ~ .x[,- 1, drop = FALSE])
  raw_errors <- map(raw_errors, ~ .x[,- 1, drop = FALSE])
  n_windows <- n_windows - 1
  former_integrated_preds <- map2(former_predictions, raw_errors, ~ mapply(function(w) prediction_integration(seeds = .x[, w], raw_errors = as.data.frame(t(.y[, 1:w]))), w = 1:n_windows, SIMPLIFY = FALSE))

  former_doxed_preds <- mapply(function(f) map(former_integrated_preds[[f]], ~ doxa_filter_plus(df[, f], .x, binary_class[f], seed)), f = 1:n_feat, SIMPLIFY = FALSE)


  windowed_df <- transpose(map(test_index, ~ df[1:.x, , drop = FALSE]))

  former_holdouts <- map(former_holdouts, ~ .x[,-1, drop = FALSE])
  windowed_df <- map(windowed_df, ~ tail(.x, -1))

  testing_errors <- mapply(function(f) pmap(list(former_holdouts[[f]], former_doxed_preds[[f]], windowed_df[[f]]), ~ custom_metrics(holdout = ..1, forecast = colMeans(..2), actuals = ..3, error_scale, error_benchmark, binary_class[f])), f = 1:n_feat, SIMPLIFY = FALSE)
  testing_errors <- as.data.frame(t(as.data.frame(map(testing_errors, ~ apply(Reduce(rbind, .x), 2, function(x) mean(x[is.finite(x)], na.rm = TRUE))))))
  rownames(testing_errors) <- NULL

  selected_preds <- map(former_doxed_preds, ~ .x[-n_windows])
  selected_holdouts <- map(former_holdouts, ~ .x[,-1])

  pred_scores <- round(as.data.frame(map2(selected_preds, selected_holdouts, ~ rowMeans(t(Reduce(rbind, map2(.x, .y, ~ prediction_score(.x, .y))))))), 4)
  rownames(pred_scores) <- NULL

  last_prediction <- engine(df, seq_len, k, method, kernel, deriv, mode, error_measurement = FALSE, error_scale, error_benchmark)$prediction

  integrated_pred <- map2(last_prediction, raw_errors, ~ prediction_integration(seeds = .x, raw_errors = as.data.frame(t(.y))))
  doxed_pred <- pmap(list(df, integrated_pred, binary_class), ~ doxa_filter_plus(..1, ..2, ..3, seed))
  prediction <- pmap(list(df, doxed_pred, binary_class), ~ fast_qpred(raw_pred = ..2, ts = ..1, ci, error_scale, error_benchmark, binary_class = ..3, dates, seed))
  prediction <- map2(prediction, pred_scores, ~ cbind(.x, pred_scores = .y))

  plot <- pmap(list(df, prediction, feat_names), ~ plotter(quant_pred = ..2, ci, ts = ..1, dates, feat_name = ..3))

  outcome <- list(prediction = prediction, testing_errors = testing_errors, plot = plot)

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
engine <- function(df, seq_len, k, method, kernel, error_measurement, deriv, mode, error_scale, error_benchmark)
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

###
fast_qpred <- function(raw_pred, ts, ci, error_scale = "naive", error_benchmark = "naive", binary_class = F, dates, seed = 42)
{
  set.seed(seed)

  quants <- sort(unique(c((1-ci)/2, 0.25, 0.5, 0.75, ci+(1-ci)/2)))

  if(binary_class == F)
  {
    p_stats <- function(x){c(min = suppressWarnings(min(x, na.rm = TRUE)), quantile(x, probs = quants, na.rm = TRUE), max = suppressWarnings(max(x, na.rm = TRUE)), mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE), mode = suppressWarnings(mlv1(x[is.finite(x)], method = "shorth")), kurtosis = suppressWarnings(kurtosis(x[is.finite(x)], na.rm = TRUE)), skewness = suppressWarnings(skewness(x[is.finite(x)], na.rm = TRUE)))}
    quant_pred <- as.data.frame(t(as.data.frame(apply(raw_pred, 2, p_stats))))
    p_value <- apply(raw_pred, 2, function(x) ecdf(x)(seq(min(raw_pred), max(raw_pred), length.out = 1000)))
    divergence <- c(max(p_value[,1] - seq(0, 1, length.out = 1000)), apply(p_value[,-1, drop = FALSE] - p_value[,-ncol(p_value), drop = FALSE], 2, function(x) abs(max(x, na.rm = TRUE))))
    upside_prob <- c(mean((raw_pred[,1]/tail(ts, 1)) > 1, na.rm = T), apply(apply(raw_pred[,-1, drop = FALSE]/raw_pred[,-ncol(raw_pred), drop = FALSE], 2, function(x) x > 1), 2, mean, na.rm = T))
    iqr_to_range <- (quant_pred[, "75%"] - quant_pred[, "25%"])/(quant_pred[, "max"] - quant_pred[, "min"])
    above_to_below_range <- (quant_pred[, "max"] - quant_pred[, "50%"])/(quant_pred[, "50%"] - quant_pred[, "min"])
    quant_pred <- round(cbind(quant_pred, iqr_to_range, above_to_below_range, upside_prob, divergence), 4)
  }

  if(binary_class == T)
  {
    p_stats <- function(x){c(min = suppressWarnings(min(x, na.rm = TRUE)), quantile(x, probs = quants, na.rm = TRUE), max = suppressWarnings(max(x, na.rm = TRUE)), prop = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE), entropy = entropy(x))}
    quant_pred <- as.data.frame(t(as.data.frame(apply(raw_pred, 2, p_stats))))
    p_value <- apply(raw_pred, 2, function(x) ecdf(x)(c(0, 1)))
    divergence <- c(max(p_value[,1] - c(0, 1)), apply(p_value[,-1, drop = FALSE] - p_value[,-ncol(p_value), drop = FALSE], 2, function(x) abs(max(x, na.rm = TRUE))))
    upgrade_prob <- c(mean(((raw_pred[,1] + 1)/tail(ts + 1, 1)) > 1, na.rm = T), apply(apply((raw_pred[,-1, drop = FALSE] + 1)/(raw_pred[,-ncol(raw_pred), drop = FALSE] + 1), 2, function(x) x > 1), 2, mean, na.rm = T))
    quant_pred <- round(cbind(quant_pred, upgrade_prob = upgrade_prob, divergence = divergence), 4)
  }


  if(is.Date(dates))
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
plotter <- function(quant_pred, ci, ts, dates = NULL, feat_name)
{
  seq_len <- nrow(quant_pred)
  n_ts <- length(ts)

  if(is.Date(dates))
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
  x_lab <- paste0("Forecasting Horizon for sequence n = ", seq_len)
  y_lab <- paste0("Forecasting Values for ", feat_name)

  lower_b <- paste0((1-ci)/2 * 100, "%")
  upper_b <- paste0((ci+(1-ci)/2) * 100, "%")

  plot <- ts_graph(x_hist = x_hist, y_hist = ts, x_forcat = x_forcat, y_forcat = quant_pred[, "50%"], lower = quant_pred[, lower_b], upper = quant_pred[, upper_b], label_x = x_lab, label_y = y_lab)
  return(plot)
}

###
doxa_filter_plus <- function(ts, mat, binary_class = F, seed)
{
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

  if(all_positive_check){mat[mat < 0] <- 0}
  if(all_negative_check){mat[mat > 0] <- 0}
  if(discrete_check){mat <- round(mat)}
  if(monotonic_increase_check){mat <- t(apply(mat, 1, function(x) monotonic_fixer(x, mode = 0)))}
  if(monotonic_decrease_check){mat <- t(apply(mat, 1, function(x) monotonic_fixer(x, mode = 1)))}

  if(binary_class == T){mat[mat > 1] <- 1; mat[mat < 1] <- 0}

  return(mat)
}


###
custom_metrics <- function(holdout, forecast, actuals, error_scale = "naive", error_benchmark = "naive", binary_class = F)
{

  if(binary_class == F)
  {
    scale <- switch(error_scale, "deviation" = sd(actuals), "naive" = mean(abs(diff(actuals))))
    benchmark <- switch(error_benchmark, "average" = rep(mean(forecast), length(forecast)), "naive" = rep(tail(actuals, 1), length(forecast)))
    me <- ME(holdout, forecast, na.rm = TRUE)
    mae <- MAE(holdout, forecast, na.rm = TRUE)
    mse <- MSE(holdout, forecast, na.rm = TRUE)
    rmsse <- RMSSE(holdout, forecast, scale, na.rm = TRUE)
    mre <- MRE(holdout, forecast, na.rm = TRUE)
    mpe <- MPE(holdout, forecast, na.rm = TRUE)
    mape <- MAPE(holdout, forecast, na.rm = TRUE)
    rmae <- rMAE(holdout, forecast, benchmark, na.rm = TRUE)
    rrmse <- rRMSE(holdout, forecast, benchmark, na.rm = TRUE)
    rame <- rAME(holdout, forecast, benchmark, na.rm = TRUE)
    mase <- MASE(holdout, forecast, scale, na.rm = TRUE)
    smse <- sMSE(holdout, forecast, scale, na.rm = TRUE)
    sce <- sCE(holdout, forecast, scale, na.rm = TRUE)
    gmrae <- GMRAE(holdout, forecast, benchmark, na.rm = TRUE)
    out <- round(c(me = me, mae = mae, mse = mse, rmsse = rmsse, mpe = mpe, mape = mape, rmae = rmae, rrmse = rrmse, rame = rame, mase = mase, smse = smse, sce = sce, gmrae = gmrae), 3)
  }

  if(binary_class == T)
  {
    dice <- suppressMessages(distance(rbind(holdout, forecast), method = "dice"))
    jaccard <- suppressMessages(distance(rbind(holdout, forecast), method = "jaccard"))
    cosine <- suppressMessages(distance(rbind(holdout, forecast), method = "cosine"))
    canberra <- suppressMessages(distance(rbind(holdout, forecast), method = "canberra"))
    gower <- suppressMessages(distance(rbind(holdout, forecast), method = "gower"))
    tanimoto <- suppressMessages(distance(rbind(holdout, forecast), method = "tanimoto"))
    hassebrook <- 1 - suppressMessages(distance(rbind(holdout, forecast), method = "hassebrook"))
    taneja <- suppressMessages(distance(rbind(holdout, forecast), method = "taneja"))
    lorentzian <- suppressMessages(distance(rbind(holdout, forecast), method = "lorentzian"))
    clark <- suppressMessages(distance(rbind(holdout, forecast), method = "clark"))
    sorensen <- suppressMessages(distance(rbind(holdout, forecast), method = "sorensen"))
    harmonic_mean <- suppressMessages(distance(rbind(holdout, forecast), method = "harmonic_mean"))
    avg <- suppressMessages(distance(rbind(holdout, forecast), method = "avg"))

    out <- round(c(dice, jaccard, cosine, canberra, gower, tanimoto, hassebrook, taneja, lorentzian, clark, sorensen, harmonic_mean, avg), 4)
  }

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

  all_data <- data.frame(x_all = c(x_hist, x_forcat), y_all = as.numeric(c(y_hist, y_forcat)))
  forcat_data <- data.frame(x_forcat = x_forcat, y_forcat = as.numeric(y_forcat))

  if(!is.null(lower) & !is.null(upper)){forcat_data$lower <- as.numeric(lower); forcat_data$upper <- as.numeric(upper)}

  plot <- ggplot()+ geom_line(data = all_data, aes_string(x = "x_all", y = "y_all"), color = hist_line, size = line_size)
  if(!is.null(lower) & !is.null(upper)){plot <- plot + geom_ribbon(data = forcat_data, aes_string(x = "x_forcat", ymin = "lower", ymax = "upper"), alpha = 0.3, fill = forcat_band)}
  plot <- plot + geom_line(data = forcat_data, aes_string(x = "x_forcat", y = "y_forcat"), color = forcat_line, size = line_size)
  if(!is.null(dbreak)){plot <- plot + scale_x_date(name = paste0("\n", label_x), date_breaks = dbreak, date_labels = date_format)}
  if(is.null(dbreak)){plot <- plot + xlab(label_x)}
  if(is.null(n_class)){plot <- plot + scale_y_continuous(name = paste0(label_y, "\n"), labels = number)}
  if(is.numeric(n_class)){plot <- plot + scale_y_continuous(name = paste0(label_y, "\n"), breaks = 1:n_class, labels = class_levels)}
  plot <- plot + ylab(label_y)  + theme_bw()
  plot <- plot + theme(axis.text=element_text(size=label_size), axis.title=element_text(size=label_size + 2))

  return(plot)
}

###
prediction_integration <- function(seeds, raw_errors){as.matrix(as.data.frame(map2(seeds, raw_errors, ~ .x + sample(.y, size = 1000, replace = TRUE))))}

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


