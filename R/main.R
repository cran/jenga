#' jenga
#'
#' @param ts A data frame with time features on columns
#' @param seq_len Positive integer. Time-step number of the projected sequence
#' @param k Positive integer. Number of neighbors to consider when applying kernel average. Min number is 3. Default: NULL (automatic selection).
#' @param method Positive integer. Distance method for calculating neighbors. Possibile options are: "euclidean", "manhattan", "chebyshev", "sorensen", "gower", "soergel", "kulczynski_d", "canberra", "lorentzian", "intersection", "non-intersection", "wavehedges", "czekanowski", "motyka", "tanimoto", "ruzicka", "inner_product", "harmonic_mean", "cosine", "hassebrook", "jaccard", "dice",  "fidelity",  "bhattacharyya", "squared_chord", "squared_euclidean", "pearson", "neyman", "squared_chi", "prob_symm", "divergence", "clark", "additive_symm", "taneja", "kumar-johnson", "avg". Default: NULL (automatic selection).
#' @param kernel String. DIstribution used to calculate kernel densities. Possible options are: "empirical", "norm", "cauchy", "logis", "unif", "t", "exp", "lnorm". Default: NULL (automatic selection).
#' @param ci Confidence interval. Default: 0.8
#' @param deriv Integer vector. Number of differentiation operations to perform for each original time feature. 0 = no change; 1: one diff; 2: two diff.
#' @param n_windows Positive integer. Number of validation tests to measure/sample error. Default: 10.
#' @param mode String. Sequencing method: deterministic ("segmented"), or non-deterministic ("sampled"). Default: NULL (automatic selection).
#' @param n_sample Positive integer. Number of samples for grid or random search. Default: 30.
#' @param search String. Two option available: "grid", "random". Default: "random".
#' @param dates Date. Vector with dates for time features.
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
#' \item errors: training and testing errors for one-step and sequence for each ts feature (rmse, mae, mdae, mpe, mape, smape)
#' \item predictions: min, max, q25, q50, q75, quantiles at selected ci, mean, sd for each ts feature
#' \item pred_stats: for each predicted time feature, IQR to range, Kullback-Leibler Divergence (compared to previous point in time), upside probability (compared to previous point in time), both averaged across all points in time and compared between the terminal and the first point in the prediction sequence.
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
#' @importFrom stats weighted.mean
#' @importFrom imputeTS na_kalman


#'@examples
#'jenga(covid_in_europe[, c(2, 3)], n_sample = 3)
#'jenga(covid_in_europe[, c(4, 5)], n_sample = 3)
#'

jenga <- function(ts, seq_len = NULL, k = NULL, method = NULL, kernel = NULL, ci = 0.8, deriv = NULL, n_windows = 10, mode = NULL, n_sample = 30, search = "random", dates = NULL, seed = 42)
{
  set.seed(seed)

  tic.clearlog()
  tic("time")

  if(!is.data.frame(ts)){stop("time features must be in dataframe format")}
  if(anyNA(ts)){ts <- as.data.frame(na_kalman(ts))}

  n_ts <- nrow(ts)
  n_feat <- ncol(ts)

  if(!is.null(k) && any(k < 3)){k[k < 3] <- 3; message("setting minimum k value to 3")}
  if(!is.null(seq_len) && any(seq_len < 3)){seq_len[seq_len < 3] <- 3; message("setting minimum seq_len value to 3")}
  if(!is.null(deriv) && any(deriv > 2)){deriv[deriv > 2] <- 2; message("setting maximum deriv value to 2")}

  if(is.null(seq_len)){sl_sample_field <- 3:floor(sqrt(n_ts))} else {sl_sample_field <- seq_len}
  if(is.null(k)){k_sample_field <- 3:floor(sqrt(n_ts))} else {k_sample_field <- k}
  if(is.null(deriv)){d_sample_field <- 0:2} else {d_sample_field <- deriv}

  if(is.null(method)){m_sample_field <- c("euclidean", "manhattan", "chebyshev", "sorensen", "gower", "soergel", "kulczynski_d", "canberra", "lorentzian", "intersection", "non-intersection", "wavehedges", "czekanowski",
                                          "motyka", "tanimoto", "ruzicka", "inner_product", "harmonic_mean", "cosine", "hassebrook", "jaccard", "dice",  "fidelity",  "bhattacharyya", "squared_chord", "squared_euclidean",
                                          "pearson", "neyman", "squared_chi", "prob_symm", "divergence", "clark", "additive_symm", "taneja", "kumar-johnson", "avg")} else {m_sample_field <- method}

  if(is.null(kernel)){krnl_sample_field <- c("norm", "cauchy", "logis", "unif", "t", "exp", "lnorm")} else {krnl_sample_field <- kernel}

  if(is.null(mode)){mode_sample_field <- c("sampled", "segmented")} else {mode_sample_field <- mode}

  n_comb <- length(sl_sample_field) * length(k_sample_field) * length(d_sample_field)^n_feat * length(m_sample_field) * length(krnl_sample_field) * length(mode_sample_field)

  if(search == "random")
  {
    if(length(sl_sample_field) == 1){sl_sample <- rep(sl_sample_field, n_sample)} else {sl_sample <- sample(sl_sample_field, n_sample, replace = TRUE)}
    if(length(k_sample_field) == 1){k_sample <- rep(k_sample_field, n_sample)} else {k_sample <- sample(k_sample_field, n_sample, replace = TRUE)}
    if(length(d_sample_field) == 1){d_sample <- replicate(n_sample, rep(d_sample_field, n_feat), simplify = F)} else {d_sample <- replicate(n_sample, sample(d_sample_field, n_feat, replace = TRUE), simplify = F)}
    if(length(d_sample_field) == n_feat){d_sample <- replicate(n_sample, d_sample_field, simplify = F)}
    if(length(m_sample_field) == 1){m_sample <- replicate(n_sample, rep(m_sample_field, n_feat), simplify = F)} else {m_sample <- replicate(n_sample, sample(m_sample_field, n_feat, replace = TRUE), simplify = F)}
    if(length(krnl_sample_field) == 1){krnl_sample <- replicate(n_sample, rep(krnl_sample_field, n_feat), simplify = F)} else {krnl_sample <- replicate(n_sample, sample(krnl_sample_field, n_feat, replace = TRUE), simplify = F)}
    if(length(mode_sample_field) == 1){mode_sample <- sample(mode_sample_field, n_sample, replace = TRUE)} else {mode_sample <- sample(mode_sample_field, n_sample, replace = TRUE)}
  }

  if(search == "grid" | n_sample > n_comb)
  {
    d_perm <- narray::split(expand.grid(replicate(n_feat, d_sample_field, simplify = F)), along = 1)
    m_perm <- narray::split(expand.grid(replicate(n_feat, m_sample_field, simplify = F)), along = 1)
    krnl_perm <- narray::split(expand.grid(replicate(n_feat, krnl_sample_field, simplify = F)), along = 1)
    grid_list <- narray::split(expand.grid(sl_sample_field, k_sample_field, d_perm, m_perm, krnl_perm, mode_sample_field), 2)
    sl_sample <- grid_list[[1]]
    k_sample <- grid_list[[2]]
    d_sample <- map(grid_list[[3]], ~ unlist(.x))
    m_sample <- map(grid_list[[4]], ~ unlist(.x))
    krnl_sample <- map(grid_list[[5]], ~ unlist(.x))
    mode_sample <- grid_list[[6]]
  }

  sl_sample <- pmap(list(sl_sample, k_sample, d_sample), ~ ((..1 * ..2) > (n_ts - max(..3) - n_windows)) * floor((n_ts - max(..3) - n_windows)/..2) +  ((..1 * ..2) <= (n_ts - max(..3) - n_windows)) * ..1)
  sl_sample <- map2(sl_sample, d_sample, ~ ((.x - max(.y)) == 0) * (.x + 2) + ((.x - max(.y)) == 1) * (.x + 1) + ((.x - max(.y)) > 1) * .x)

  if(length(seq_len)==1 & any(sl_sample != seq_len)){message("fixing seq_len for available data")}

  hyper_params <- list(sl_sample, k_sample, d_sample, m_sample, krnl_sample, mode_sample)

  exploration <- suppressMessages(pmap(hyper_params, ~ tryCatch(hood(ts, seq_len = ..1, k = ..2, method = ..4, kernel = ..5, ci, deriv = ..3, n_windows, mode = ..6, dates = dates), error = function(e) NA)))
  not_na <- !is.na(exploration)
  #failures <- Reduce(cbind, list(seq_len = sl_sample[!not_na], k = k_sample[!not_na], deriv = d_sample[!not_na], dist_method = m_sample[!not_na], kernel = krnl_sample[!not_na], mode = mode_sample[!not_na]))
  exploration <- exploration[not_na]

  if(n_feat > 1){test_metrics <- map_dfr(exploration, ~ apply(Reduce(rbind, .x$test_metrics), 2, mean))}
  if(n_feat == 1){test_metrics <- map_dfr(exploration, ~ unlist(.x$test_metrics))}
  weights <- apply(test_metrics, 2, function(x) {sd(x[is.finite(x)], na.rm = TRUE)/abs(mean(x[is.finite(x)], na.rm = TRUE))})
  finite_w <- is.finite(weights)

  history <- Reduce(cbind, list(seq_len = sl_sample[not_na], k = k_sample[not_na], deriv = d_sample[not_na], dist_method = m_sample[not_na], kernel = krnl_sample[not_na], mode = mode_sample[not_na]))
  colnames(history) <- c("seq_len", "k", "deriv", "dist_method", "kernel", "mode")
  history <- cbind(history, test_metrics)
  history$wgt_avg_rank <- round(apply(apply(test_metrics[, finite_w, drop = FALSE], 2, rank), 1, weighted.mean, w = weights[finite_w]), 2)
  history <- history[order(history$wgt_avg_rank),]

  best_index <- as.numeric(rownames(history[1,]))

  best <- list(wt_avg_best = history[1, c(1:6, 15)], model = exploration[[best_index]])

  toc(log = TRUE)
  time_log <- seconds_to_period(round(parse_number(unlist(tic.log())), 0))

  outcome <- list(exploration = exploration, history = history, best = best, time_log = time_log)

  return(outcome)
}

