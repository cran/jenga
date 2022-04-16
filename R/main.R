#' jenga
#'
#' @param ts A data frame with time features on columns
#' @param seq_len Positive integer. Time-step number of the projected sequence
#' @param k Positive integer. Number of neighbors to consider when applying kernel average. Min number is 3. Default: NULL (automatic selection).
#' @param method Positive integer. Distance method for calculating neighbors. Possibile options are: "euclidean", "manhattan", "chebyshev", "gower", "lorentzian", "jaccard", "dice", "squared_euclidean", "divergence", "clark", "avg". Default: NULL (automatic selection).
#' @param kernel String. Distribution used to calculate kernel densities. Possible options are: "norm", "cauchy", "logis", "unif", "t", "lnorm". Default: NULL (automatic selection).
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
#' \item best: results for the best model in term of weighted average rank, including:
#' \itemize{
#' \item errors: training and testing errors for one-step and sequence for each ts feature (me, mae, mse, rmsse, mpe, mape, rmae, rrmse, rame, mase, smse, sce, gmrae, pred_score)
#' \item predictions: min, max, q25, q50, q75, quantiles at selected ci, mean, sd, mode, skewness, kurtosis, IQR to range, risk ratio, upside probability, divergence and entropy for each point fo predicted sequences
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
#' @importFrom stats weighted.mean lm
#' @importFrom imputeTS na_kalman


#'@examples
#'jenga(covid_in_europe[, c(2, 3)], n_sample = 1)
#'jenga(covid_in_europe[, c(4, 5)], n_sample = 1)
#'

jenga <- function(ts, seq_len = NULL, k = NULL, method = NULL, kernel = NULL, ci = 0.8, n_windows = 10, mode = NULL, n_sample = 30, search = "random", dates = NULL, error_scale = "naive", error_benchmark = "naive", seed = 42)
{
  ###SUPPORT
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
  set.seed(seed)

  tic.clearlog()
  tic("time")

  if(!is.data.frame(ts)){stop("time features must be in dataframe format")}
  if(anyNA(ts)){ts <- as.data.frame(na_kalman(ts))}

  n_ts <- nrow(ts)
  n_feat <- ncol(ts)
  deriv <- map_dbl(ts, ~ best_deriv(.x))
  fixed <- FALSE
  if(all(c(length(seq_len)==1, length(k)==1, length(method)<=n_feat, length(kernel)<=n_feat, length(mode)==1))){fixed <- TRUE}

  if(!is.null(k) && any(k < 3)){k[k < 3] <- 3; message("setting minimum k value to 3")}
  if(!is.null(seq_len) && any(seq_len < 3)){seq_len[seq_len < 3] <- 3; message("setting minimum seq_len value to 3")}

  ###if(fixed == TRUE & (is.null(seq_len) | is.null(k) |is.null(method) | is.null(kernel))){stop("fixed flag requires values for seq_len, k, method and kernel")}

  if(fixed == FALSE)
  {
  if(is.null(seq_len)){sl_sample_field <- 3:floor(sqrt(n_ts))} else {sl_sample_field <- seq_len}
  if(is.null(k)){k_sample_field <- 3:floor(sqrt(n_ts))} else {k_sample_field <- k}

  if(is.null(method)){m_sample_field <- c("euclidean", "manhattan", "chebyshev", "gower", "lorentzian", "jaccard", "dice", "squared_euclidean", "divergence", "clark", "avg")} else {m_sample_field <- method}

  ####OLD "euclidean", "manhattan", "chebyshev", "sorensen", "gower", "soergel", "kulczynski_d", "canberra", "lorentzian", "wavehedges", "czekanowski", "tanimoto", "jaccard", "dice",  "squared_chord", "squared_euclidean", "pearson", "neyman", "squared_chi", "prob_symm", "divergence", "clark", "additive_symm", "taneja", "kumar-johnson", "avg"

  if(is.null(kernel)){krnl_sample_field <- c("norm", "cauchy", "logis", "unif", "t")} else {krnl_sample_field <- kernel}

  if(is.null(mode)){mode_sample_field <- c("sampled", "segmented")} else {mode_sample_field <- mode}

  n_comb <- length(sl_sample_field) * length(k_sample_field) * length(m_sample_field)^n_feat * length(krnl_sample_field)^n_feat * length(mode_sample_field)

  if(search == "random")
  {
    if(length(sl_sample_field) == 1){sl_sample <- rep(sl_sample_field, n_sample)} else {sl_sample <- sample(sl_sample_field, n_sample, replace = TRUE)}
    if(length(k_sample_field) == 1){k_sample <- rep(k_sample_field, n_sample)} else {k_sample <- sample(k_sample_field, n_sample, replace = TRUE)}
    if(length(m_sample_field) == 1){m_sample <- replicate(n_sample, rep(m_sample_field, n_feat), simplify = F)} else {m_sample <- replicate(n_sample, sample(m_sample_field, n_feat, replace = TRUE), simplify = F)}
    if(length(krnl_sample_field) == 1){krnl_sample <- replicate(n_sample, rep(krnl_sample_field, n_feat), simplify = F)} else {krnl_sample <- replicate(n_sample, sample(krnl_sample_field, n_feat, replace = TRUE), simplify = F)}
    if(length(mode_sample_field) == 1){mode_sample <- sample(mode_sample_field, n_sample, replace = TRUE)} else {mode_sample <- sample(mode_sample_field, n_sample, replace = TRUE)}
  }

  if(search == "grid" | n_sample > n_comb)
  {
    m_perm <- narray::split(expand.grid(replicate(n_feat, m_sample_field, simplify = F)), along = 1)
    krnl_perm <- narray::split(expand.grid(replicate(n_feat, krnl_sample_field, simplify = F)), along = 1)
    grid_list <- narray::split(expand.grid(sl_sample_field, k_sample_field, m_perm, krnl_perm, mode_sample_field), 2)
    sl_sample <- grid_list[[1]]
    k_sample <- grid_list[[2]]
    m_sample <- map(grid_list[[3]], ~ unlist(.x))
    krnl_sample <- map(grid_list[[4]], ~ unlist(.x))
    mode_sample <- grid_list[[5]]
  }

  sl_sample <- pmap_dbl(list(sl_sample, k_sample, replicate(n_sample, deriv, simplify = F)), ~ ((..1 * ..2) > (n_ts - max(..3) - n_windows)) * floor((n_ts - max(..3) - n_windows)/..2) +  ((..1 * ..2) <= (n_ts - max(..3) - n_windows)) * ..1)
  sl_sample <- map2_dbl(sl_sample, replicate(n_sample, deriv, simplify = F), ~ ((.x - max(.y)) == 0) * (.x + 2) + ((.x - max(.y)) == 1) * (.x + 1) + ((.x - max(.y)) > 1) * .x)

  if(length(seq_len)==1 && any(sl_sample != seq_len)){message("fixing seq_len for available data")}

  hyper_params <- list(sl_sample, k_sample, m_sample, krnl_sample, mode_sample)
  }

  if(fixed == TRUE){hyper_params <- list(sl_sample = list(seq_len), k_sample = list(k), m_sample = list(method),  krnl_sample = list(kernel), mode_sample = list(mode))}

  exploration <- suppressMessages(pmap(hyper_params, ~ tryCatch(suppressWarnings(hood(ts, seq_len = ..1, k = ..2, method = ..3, kernel = ..4, ci, deriv = deriv, n_windows, mode = ..5, dates = dates, error_scale, error_benchmark)), error = function(e) NA)))
  not_na <- !is.na(exploration)
  exploration <- exploration[not_na]

  #return(map(hyper_params, ~ .x[!not_na]))


  if(fixed == FALSE)
  {
  test_metrics <- map_dfr(exploration, ~ apply(.x$test_metrics, 2, function(x) mean(x)))
  colnames(test_metrics) <- paste0("avg_", colnames(test_metrics))
  pred_scores <- map_dbl(exploration, ~ min(unlist(.x$test_metrics$pred_score)))
  history <- Reduce(cbind, list(seq_len = sl_sample[not_na], k = k_sample[not_na], dist_method = m_sample[not_na], kernel = krnl_sample[not_na], mode = mode_sample[not_na]))
  colnames(history) <- c("seq_len", "k", "dist_method", "kernel", "mode")
  history <- as.data.frame(cbind(history, pred_scores, test_metrics))
  if(n_sample > 1){history <- history[order(-unlist(history$pred_scores), abs(unlist(history$avg_me)), unlist(history$avg_mae)),]}
  best_index <- as.numeric(rownames(history[1,]))

  best <- exploration[[best_index]]
  }

  if(fixed == TRUE)
  {
  history <- NULL
  best <- exploration[[1]]
  }

  toc(log = TRUE)
  time_log <- seconds_to_period(round(parse_number(unlist(tic.log())), 0))

  outcome <- list(exploration = exploration, history = history, best = best, time_log = time_log)

  return(outcome)
}

