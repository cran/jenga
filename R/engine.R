#' support functions for jenga
#'
#' @author Giancarlo Vercellino \email{giancarlo.vercellino@gmail.com}
#'
#' @importFrom narray split
#' @importFrom utils head tail
#' @importFrom modeest mlv1
#' @importFrom moments skewness kurtosis
#' @importFrom stats dcauchy density dexp dlnorm dlogis dnorm dt dunif median sd
#' @import purrr
#' @import abind
#' @import philentropy
#' @import abind
#' @import ggplot2
#' @import greybox

#' @param ts A data frame.
#' @param seq_len Positive integer.
#' @param k Positive integer.
#' @param method Positive integer.
#' @param kernel String.
#' @param error_measurement Logical.
#' @param deriv Integer vector.
#' @param mode String.
#' @param error_scale String.
#' @param error_benchmark String.


###
engine <- function(ts, seq_len, k, method, kernel, error_measurement, deriv, mode, error_scale, error_benchmark)
{

  my_metrics <- function(holdout, forecast, past, error_scale, error_benchmark)
  {
    scale <- switch(error_scale, "deviation" = sd(past), "naive" = mean(abs(diff(past))))
    benchmark <- switch(error_benchmark, "average" = rep(mean(forecast), length(forecast)), "naive" = rep(tail(past, 1), length(forecast)))

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
    return(out)
  }

  segmenter <- function(vector, seq_len)
  {
    if(seq_len > length(vector)){stop("vector too short for sequence length")}
    n_segment <- floor(length(vector)/seq_len)
    residual <- length(vector)%%seq_len
    fixed_vector <- vector
    if(residual > 0){fixed_vector <- tail(vector, - residual)}
    segment_id <- sort(rep(1:n_segment, seq_len))
    seq_list <- narray::split(fixed_vector, along = 1, subsets = segment_id)
    segmented <- as.data.frame(Reduce(rbind, seq_list))
    rownames(segmented) <- NULL
    return(segmented)
  }

  reframe<-function(data, length)
  {
    if(!is.data.frame(data)){data <- data.frame(data)}
    slice_list <- narray::split(data, along=2)
    reframed <- abind(map(slice_list, ~ t(apply(embed(.x, dimension=length), 1, rev))), along=3)
    return(reframed)
  }

  recursive_diff <- function(vector, deriv)
  {
    head_value <- vector("numeric", deriv)
    tail_value <- vector("numeric", deriv)
    if(deriv==0){head_value = NULL; tail_value = NULL}
    if(deriv > 0){for(i in 1:deriv){head_value[i] <- head(vector, 1); tail_value[i] <- tail(vector, 1); vector <- diff(vector)}}
    outcome <- list(vector = vector, head_value = head_value, tail_value = tail_value)
    return(outcome)
  }

  invdiff <- function(vector, heads)
  {
    if(is.null(heads)){return(vector)}
    for(d in length(heads):1){vector <- cumsum(c(heads[d], vector))}
    return(vector)
  }

  orig <- ts
  n_feat <- ncol(orig)
  feat_names <- colnames(orig)

  diff_model <- map2(orig, deriv, ~ recursive_diff(.x, .y))
  ts_list <- map(diff_model, ~ .x$vector)

  if(mode == "segmented")
  {
    sequenced_list <- map2(ts_list, deriv, ~ segmenter(.x, seq_len - .y))
    min_size <- min(map_dbl(sequenced_list, ~ nrow(.x)))
    sequenced_list <- map(sequenced_list, ~ tail(.x, min_size))
  }

  if(mode == "sampled")
  {
    fixed_ts_list <- map2(ts_list, deriv, ~ head(.x, - (seq_len - .y)))
    sequenced_list <- map2(fixed_ts_list, deriv, ~ as.data.frame(reframe(.x, seq_len - .y)))
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

  distances <- map2(sequenced_list, method, ~ mapply(function(x) suppressMessages(distance(rbind(.x[n_seq,], .x[x,]), method = .y)), x = 1:n_seq))

  norm_distances <- map_dfc(distances, ~ .x/sum(.x))
  dmatrix_line <- apply(norm_distances, 1, mean)

  ranking <- rank(dmatrix_line)
  seq_id <- 1:n_seq
  neighbour_id <- seq_id[ranking <= k]
  #neighbour_id <- seq_id[ranking <= (k + 1) & !ranking==1] + 1
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
  prediction <- suppressMessages(map2_dfc(prediction_list, diff_model, ~ tail(invdiff(.x, heads = .y$tail_value), seq_len)))
  colnames(prediction) <- feat_names

  raw_error <- NULL
  test_metrics <- NULL
  if(error_measurement == TRUE)
  {
    raw_error <- actual - prediction
    test_metrics <- pmap(list(actual, prediction, past), ~ my_metrics(..1, ..2, ..3, error_scale, error_benchmark))
  }

  outcome <- list(prediction = prediction, actual = actual, raw_error = raw_error, test_metrics = test_metrics)
  return(outcome)
}
