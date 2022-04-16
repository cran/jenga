#' support functions for jenga
#'
#' @author Giancarlo Vercellino \email{giancarlo.vercellino@gmail.com}
#'
#' @importFrom narray split
#' @importFrom modeest mlv1
#' @importFrom moments skewness kurtosis
#' @importFrom stats quantile sd ecdf
#' @importFrom scales number
#' @importFrom lubridate is.Date
#' @import purrr
#' @import entropy

#' @param ts A data frame.
#' @param seq_len Positive integer.
#' @param k Positive integer.
#' @param method Positive integer.
#' @param kernel String.
#' @param ci Confidence interval.
#' @param deriv Integer vector.
#' @param n_windows Positive integer.
#' @param mode String.
#' @param dates Date.
#' @param error_scale String.
#' @param error_benchmark String.

  hood <- function(ts, seq_len, k, method, kernel, ci, deriv, n_windows, mode, dates, error_scale, error_benchmark)
  {
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

    prediction_score <- function(integrated_preds, ground_truth)
    {
      pfuns <- apply(integrated_preds, 2, ecdf)
      pvalues <- purrr::map2_dbl(pfuns, ground_truth, ~ .x(.y))
      scores <- mean(1 - 2 * abs(pvalues - 0.5))
      return(scores)
    }

    doxa_filter <- function(orig, mat, n_class = NULL)
    {
      if(is.list(orig)){orig <- unlist(orig)}
      if(!is.matrix(mat)){mat <- as.matrix(mat)}
      discrete_check <- all(orig%%1 == 0)
      all_positive_check <- all(orig >= 0)
      all_negative_check <- all(orig <= 0)
      monotonic_increase_check <- all(diff(orig) >= 0)
      monotonic_decrease_check <- all(diff(orig) <= 0)
      class_check <- FALSE
      if(is.integer(n_class)){class_check <- length(unique(orig)) <= n_class}

      monotonic_fixer <- function(x, mode)
      {
        model <- recursive_diff(x, 1)
        vect <- model$vector
        if(mode == 0){vect[vect < 0] <- 0; vect <- invdiff(vect, model$head_value, add = T)}
        if(mode == 1){vect[vect > 0] <- 0; vect <- invdiff(vect, model$head_value, add = T)}
        return(vect)
      }

      if(discrete_check){mat <- floor(mat)}
      if(monotonic_increase_check){mat <- t(apply(mat, 1, function(x) monotonic_fixer(x, mode = 0)))}
      if(monotonic_decrease_check){mat <- t(apply(mat, 1, function(x) monotonic_fixer(x, mode = 1)))}
      if(all_positive_check){mat[mat < 0] <- 0}
      if(all_negative_check){mat[mat > 0] <- 0}
      if(class_check){mat[!(mat %in% unique(orig))] <- ((mat[!(mat %in% unique(orig))] > max(unique(orig))) * max(unique(orig))) + ((mat[!(mat %in% unique(orig))] < min(unique(orig))) * min(unique(orig)))}

      df <- as.data.frame(mat)
      return(df)
    }

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

    invdiff <- function(vector, heads, add = FALSE)
    {
      vector <- unlist(vector)
      if(is.null(heads)){return(vector)}
      for(d in length(heads):1){vector <- cumsum(c(heads[d], vector))}
      if(add == FALSE){return(vector[-c(1:length(heads))])} else {return(vector)}
    }

    ###
    ts <- as.data.frame(ts)
    n_feat <- ncol(ts)
    n_ts <- nrow(ts)
    feat_names <- colnames(ts)

    p_stats <- function(x){stats <- round(c(min = suppressWarnings(min(x, na.rm = TRUE)), quantile(x, probs = quants, na.rm = TRUE), max = suppressWarnings(max(x, na.rm = TRUE)), mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE), mode = tryCatch(suppressWarnings(mlv1(x[is.finite(x)], method = "shorth")), error = function(e) NA), kurtosis = tryCatch(suppressWarnings(kurtosis(x[is.finite(x)], na.rm = TRUE)), error = function(e) NA), skewness = tryCatch(suppressWarnings(skewness(x[is.finite(x)], na.rm = TRUE)), error = function(e) NA)), 3); return(stats)}

    all_positive_check <- map_lgl(ts, ~ ifelse(all(.x >= 0), TRUE, FALSE))
    all_negative_check <- map_lgl(ts, ~ ifelse(all(.x < 0), TRUE, FALSE))

    if(seq_len * k > n_ts){stop("vector length too short for testing with seq_len and k")}
    test_index <- unique(round(seq(seq_len * k, n_ts, length.out = n_windows)))
    if(length(test_index) < n_windows){message("testing on ", length(test_index), " windows")}

    results <- map(test_index, ~ tryCatch(engine(ts[1:.x,, drop = FALSE], seq_len, k, method, kernel, deriv, mode, error_measurement = TRUE, error_scale, error_benchmark), error = function(e) NA))
    not_na <- !is.na(results)
    results <- results[not_na]
    raw_errors <- map(transpose(map(results, ~ narray::split(.x$raw_error, along = 2))), ~ as.data.frame(.x))
    former_predictions <- map(transpose(map(results, ~ narray::split(.x$prediction, along = 2))), ~ as.data.frame(.x))
    former_actuals <- map(transpose(map(results, ~ narray::split(.x$actual, along = 2))), ~ as.data.frame(.x))

    pred_integration <- function(pred_seed, raw_error){as.data.frame(map2(pred_seed, raw_error, ~ .x + sample(.y, size = 1000, replace = TRUE)))}
    former_sample_pred <- map2(former_predictions, raw_errors, ~ mapply(function(i) pred_integration(.x[, i], as.data.frame(t(.y))), i = 1:ncol(.x), SIMPLIFY = FALSE))
    former_sample_pred <- lapply(1:n_feat, function(f) lapply(1:n_windows, function(w) doxa_filter(ts[, f], former_sample_pred[[f]][[w]])))
    ##former_sample_pred <- map_depth(former_sample_pred, 2, ~ as.data.frame(.x))
    names(former_sample_pred) <- feat_names

    #former_sample_pred <- map2(ts, former_sample_pred, ~ mapply(function(i) doxa_filter(.x, i), i = .y, SIMPLIFY = FALSE))

    pred_scores <- mapply(function(p, t) prediction_score(p, t), p = flatten(former_sample_pred), t = flatten(former_actuals), SIMPLIFY = TRUE)
    pred_scores <- split(pred_scores, along = 1, subsets = sort(base::rep(1:n_feat, each=n_windows)))
    avg_pred_scores <- map(pred_scores, ~ mean(.x, na.rm = TRUE))

    last_prediction <- engine(ts, seq_len, k, method, kernel, deriv, mode, error_measurement = FALSE, error_scale, error_benchmark)$prediction

    window_metrics <- map(transpose(map(results, ~ .x$test_metrics)), ~ Reduce(rbind, .x))
    if(n_feat > 1){test_metrics <- Reduce(rbind, map(window_metrics, ~ apply(.x, 2, function(x) mean(x[is.finite(x)], na.rm = TRUE))))}
    if(n_feat == 1){test_metrics <- t(as.data.frame(map(window_metrics, ~ apply(.x, 2, function(x) mean(x[is.finite(x)], na.rm = TRUE)))))}
    rownames(test_metrics) <- feat_names
    test_metrics <- as.data.frame(cbind(test_metrics, pred_score = round(unlist(avg_pred_scores), 4)))

    quants <- sort(unique(c((1-ci)/2, 0.25, 0.5, 0.75, ci+(1-ci)/2)))

    sample_pred <- map2(last_prediction, raw_errors, ~ map2_dfc(.x, narray::split(.y, along = 1), ~ .x + sample(unlist(.y), size = 1000, replace = TRUE)))
    sample_pred <- map2(ts, sample_pred, ~ doxa_filter(.x, .y))

    prediction <- map(sample_pred, ~ as.data.frame(t(apply(.x, 2, p_stats))))
    iqr_to_range <- map(prediction, ~ tryCatch((.x[, "75%"] - .x[, "25%"])/(.x[, "max"] - .x[, "min"]), error = function(e) NA))
    risk_ratio <- map(prediction, ~ tryCatch((.x[, "max"] - .x[, "50%"])/(.x[, "50%"] - .x[, "min"]), error = function(e) NA))
    growths <- map(prediction, ~ mapply(function(m, s) rnorm(1000, m, s), m = .x[, "mean"], s = .x[, "sd"]))
    upside_prob <- map(growths, ~ tryCatch(c(NA, colMeans(apply(.x[,-1]/.x[, - seq_len], 2, function(x) x > 1))), error = function(e) NA))
    pvalues <- map(sample_pred, ~ as.data.frame(map(.x, ~ ecdf(.x)(.x))))
    divergence <- map(pvalues, ~ tryCatch(c(NA, apply(.x[,-1] - .x[, - seq_len], 2, function(x) abs(max(x, na.rm = T)))), error = function(e) NA))
    entropy <- map(sample_pred, ~ tryCatch(apply(.x, 2, entropy), error = function(e) NA))
    prediction <- pmap(list(prediction, iqr_to_range, risk_ratio, upside_prob, divergence, entropy), ~ round(cbind(..1, iqr_to_range = ..2, risk_ratio = ..3, upside_prob = ..4, divergence = ..5, entropy = ..6), 4))
    names(prediction) <- feat_names
    prediction <- map(prediction, ~ {rownames(.x) <- NULL; return(.x)})

    x_hist <- 1:n_ts
    x_forcat <- (n_ts + 1):(n_ts + seq_len)

    if(is.Date(dates) & length(dates) == n_ts)
    {
      pred_dates <- tail(dates, seq_len) + seq_len
      prediction <- map(prediction, ~ {rownames(.x) <- pred_dates; return(.x)})
      x_hist <- dates
      x_forcat <- pred_dates
    }

    plot <- pmap(list(prediction, ts, feat_names), ~ ts_graph(x_hist = x_hist, y_hist = ..2, x_forcat = x_forcat, y_forcat = ..1[, 3], lower = ..1[, 1], upper = ..1[, 5], label_y = ..3))

    outcome <- list(prediction = prediction, test_metrics = test_metrics, plot = plot)

    return(outcome)
  }
